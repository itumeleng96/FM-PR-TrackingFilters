classdef kalmanFilter

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,k_d;
    end
    
    methods
        function obj = kalmanFilter(dt,std_acc,r_std,rdot_std,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % x_std_meas : standard deviation of the measurement in the x-direction
            
        
            obj.X= X_initial;                   % Initial State
            obj.dt = dt;                        % Update Interval
            c=299792458;
            obj.k_d = -c/94e6;                  % Wave number k=-lambda=c/f

            obj.F = [1, obj.k_d*dt;0, 1];       % State transition matrix
                    
            
            obj.H = [1,0;0,1;];                 % Measurement Function

            
            obj.Q = [(dt^4)/4,(dt^3)/2;  %Process uncertainty
                    (dt^3)/2,dt^2]*std_acc;

            obj.R = [r_std^2,0;0,rdot_std^2];   % Measurement Uncertainty
            obj.S = [2.5e3,0;0,0.0005];                  % System Uncertainty  
            obj.P = eye(size(obj.F,2));         % Uncertainty Covariance

        end
        
        function [X_pred,KF_obj1] = predict(obj)

            %Predict next state (prior)
            
            % x = Fx
            obj.X = obj.F*obj.X ;
            
            % P = FPF' + Q
            obj.P = obj.F * obj.P * obj.F.' + obj.Q;

            %Return Prior
            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %UPDATE STAGE
            
            %S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            disp(obj.S);
            %K = PH'inv(S)
            K = obj.P * obj.H.'*obj.S^(-1);

            %x = x + Ky
            obj.X = obj.X + K * (z-obj.H * obj.X);
            
            I = eye(size(obj.H,2));

            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P ;
            
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


