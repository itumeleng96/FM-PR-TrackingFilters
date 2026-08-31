%CSKF  Covariance-Scaling Kalman Filter for bistatic range-Doppler tracking.
%
%   State vector:   x = [r, rdot, f, fdot]'  (range, range-rate, Doppler, Doppler-rate)
%   Measurement:    z = [r, f]'              (range and Doppler centroid)
%   Motion model:   Nearly Constant Velocity (NCV) with sampling interval dt
%
%   Adaptive covariance scaling (Akhlaghi, Zhou & Huang, 2018) is applied to
%   both R and Q on every accepted measurement, with forgetting factor
%   obj.cs_alpha (default 0.98). Q is adapted when obj.cs_adapt_q is true.
%
%   Construction:
%     f = CSKF(dt, std_acc, r_std, rdot_std, X_initial)
%       dt          - scan interval (s)
%       std_acc     - [ax, ay] process-noise std for range and Doppler axes
%       r_std       - measurement std, range channel
%       rdot_std    - measurement std, Doppler channel
%       X_initial   - initial state vector [r0; rdot0; f0; fdot0]
%
%   Public methods:
%     [xhat, obj] = predict(obj)      one-step prediction
%     [xhat, obj] = update(obj, z)    measurement update on z = [r; f]
%
%   Adaptive tuning knobs (set after construction):
%     obj.cs_alpha    - forgetting factor in [0,1]           (default 0.98)
%     obj.cs_adapt_q  - true = also adapt Q, false = R only  (default false)
%
%   Reference:
%     Akhlaghi, S., Zhou, N., Huang, Z. (2018). "Adaptive adjustment of noise
%     covariance in Kalman filter for dynamic state estimation." IEEE PES
%     General Meeting. arXiv:1702.00884.
classdef CSKF

    properties
        dt, X, F, H, Q, Q_nominal, R, R_nominal, P, S, std_acc, wk, ...
        cs_alpha, cs_adapt_q;
    end
    
    methods
        function obj = CSKF(dt,std_acc,r_std,rdot_std,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % x_std_meas : standard deviation of the measurement in the x-direction
            
        
            obj.X= X_initial;                              % Initial State
            obj.dt = dt;                                   % Update Interval

                    
            obj.F = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];                            % Measurement Function

            


            obj.Q = [std_acc(1)*(dt^4)/4,std_acc(1)*(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     std_acc(1)*(dt^3)/2, std_acc(1)*dt^2, 0, 0;
                     0, 0, std_acc(2)*(dt^4)/4,std_acc(2)*(dt^3)/2;
                     0, 0, std_acc(2)*(dt^3)/2, std_acc(2)*dt^2];

            obj.R = [r_std,0;
                     0,rdot_std];                  % Measurement Uncertainty
            obj.R_nominal = obj.R;                 % Nominal R for fading-memory decay
            
            obj.P = [5,0,0,0;                              % Initial Error Covariance Matrix
                     0, 1, 0, 0;
                     0, 0, 2,0;
                     0, 0, 0, 1]; 
            

            obj.wk = [std_acc(1)*dt^2;std_acc(1)*dt;std_acc(2)*dt^2;std_acc(2)*dt];
            obj.S = obj.R;

            % Akhlaghi (2018) covariance-scaling defaults (settable post-construction).
            obj.cs_alpha   = 0.98;   % forgetting factor (both R and Q)
            obj.cs_adapt_q = true;   % master switch for Q adaptation
            obj.Q_nominal  = obj.Q;  % nominal Q retained as PSD floor



        end
        
        function [X_pred,KF_obj1] = predict(obj)

            %PREDICT NEXT STATE (prior)
            % x = Fx
            obj.X = obj.F*obj.X +obj.wk; 
            
            % P = FPF' + Q
            obj.P = obj.F * obj.P * obj.F.' + obj.Q;
            X_pred = obj.X;
            KF_obj1  = obj;

        end
        
        function [X_est, KF_obj2] = update(obj, z)
            % Akhlaghi (2018) verbatim continuous R + Q. Outlier rejection is
            % handled upstream at the association step (GNN gate in runner),
            % so the filter unconditionally applies the update.
            %   R <- alpha R + (1-alpha)( eps eps' + H P^- H' )       Eq. 11
            %   Q <- alpha Q + (1-alpha)( K d d' K' )   per block     Eq. 15
            alpha  = obj.cs_alpha;
            P_pred = obj.P;

            % Innovation and standard KF update.
            obj.S  = obj.H * P_pred * obj.H.' + obj.R;
            ek     = z - obj.H * obj.X;
            K      = (P_pred * obj.H.') / obj.S;
            obj.X  = obj.X + K * ek;
            I      = eye(size(obj.H, 2));
            obj.P  = (I - K * obj.H) * P_pred;

            % Akhlaghi Eq. 11 — continuous R update (residual-based).
            eps_post = z - obj.H * obj.X;
            HPH      = obj.H * P_pred * obj.H.';
            obj.R    = alpha * obj.R + (1 - alpha) * (eps_post * eps_post.' + HPH);

            % Akhlaghi Eq. 15 — continuous Q update per block (range 1:2, doppler 3:4).
            if obj.cs_adapt_q
                KddK = K * (ek * ek.') * K.';
                obj.Q(1:2,1:2) = alpha * obj.Q(1:2,1:2) + (1 - alpha) * KddK(1:2,1:2);
                obj.Q(3:4,3:4) = alpha * obj.Q(3:4,3:4) + (1 - alpha) * KddK(3:4,3:4);
                obj.Q(1:2,1:2) = max(obj.Q(1:2,1:2), obj.Q_nominal(1:2,1:2));
                obj.Q(3:4,3:4) = max(obj.Q(3:4,3:4), obj.Q_nominal(3:4,3:4));
            end

            X_est   = obj.X;
            KF_obj2 = obj;
        end
    end

 end
        
    


