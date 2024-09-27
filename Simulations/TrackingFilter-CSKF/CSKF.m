classdef CSKF

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,wk,residuals;
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

            

            obj.Q = std_acc*[(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            obj.R = [r_std,0;
                     0,rdot_std];                  % Measurement Uncertainty
            
            obj.P = [5,0,0,0;                              % Initial Error Covariance Matrix
                     0, 1, 0, 0;
                     0, 0, 2,0;
                     0, 0, 0, 1]; 
            
            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];

            obj.wk = std_acc*[dt^2;dt;dt^2;dt];
            obj.S = obj.R;

            obj.residuals =[];


        end
        
        function [X_pred,KF_obj1] = predict(obj)

            %PREDICT NEXT STATE (prior)
            % x = Fx
            obj.X = obj.F*obj.X +obj.wk; 
            
            % P = FPF' + Q
            obj.P = obj.A * obj.P * obj.A.' + obj.Q;
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            X_pred = obj.X;
            KF_obj1  = obj;

        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %UPDATE STAGE

            %Adaptive estimation of R matrix
           

            %S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Covariance Scaling Kalman Step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mk = e'.[HPH' +R]^−1.ek .
            
            % Kalman prediction step: S = H*P*H' + R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            
            % Compute residual (measurement innovation): ek = z - H*X
            ek = z - obj.H * obj.X;
            
            % Mahalanobis distance (combined Doppler and Range)
            mahalanobisDistSquared = ek' * (obj.S \ ek);
            
            % Number of samples to average
            M = 6;
            r_adapt = obj.R;
            alpha =0.6;
            % Moving window for residuals
            obj.residuals = [obj.residuals, mahalanobisDistSquared];
            
            % Check for outliers using a combined approach for Doppler and Range
            if size(obj.residuals, 2) > M
                residual_mean = mean(obj.residuals(end-M:end-1));
                residual_std = std(obj.residuals(end-M:end-1));
                
                if mahalanobisDistSquared > 1 && mahalanobisDistSquared > (residual_mean + residual_std)                    
                    r_adapt = alpha * r_adapt + (1 - alpha) * (mahalanobisDistSquared + obj.H * obj.P * obj.H');
                    obj.residuals = obj.residuals(end-M:end-1);
                end
            end
            
            % Update R with the adapted covariance matrix
            
            S1 = obj.H * obj.P * obj.H.' + r_adapt;
            %K = PH'inv(S)
            K = (obj.P * obj.H.') * S1^(-1);

            obj.X = obj.X + K * (z-obj.H * obj.X);
            
            I = eye(size(obj.H,2));
    
            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P ; 
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


