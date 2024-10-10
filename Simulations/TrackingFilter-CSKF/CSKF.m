classdef CSKF

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,wk,residuals_x,residuals_y;
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
            
            obj.P = [5,0,0,0;                              % Initial Error Covariance Matrix
                     0, 1, 0, 0;
                     0, 0, 2,0;
                     0, 0, 0, 1]; 
            

            obj.wk = [std_acc(1)*dt^2;std_acc(1)*dt;std_acc(2)*dt^2;std_acc(2)*dt];
            obj.S = obj.R;

            obj.residuals_x =[];
            obj.residuals_y =[];



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
            %UPDATE STAGE
        
            % Adaptive estimation of R matrix
        
            % S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Covariance Scaling Kalman Step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute residual (measurement innovation): ek = z - H*X
            ek = z - obj.H * obj.X;
        
            % Separate residuals for x and y components
            ek_x = ek(1); % Residual for x
            ek_y = ek(2); % Residual for y
        
            % Mahalanobis distance (separate for x and y)
            mahalanobisDistSquared_x = ek_x' * (obj.S(1,1) \ ek_x);
            mahalanobisDistSquared_y = ek_y' * (obj.S(2,2) \ ek_y);
        
            % Number of samples to average
            M = 6;
            r_adapt_x = obj.R(1,1);  % Adaptive R for x
            r_adapt_y = obj.R(2,2);  % Adaptive R for y
            alpha = 0.7;  % Smoothing factor
        
            % Moving window for residuals (for x and y)
            obj.residuals_x = [obj.residuals_x, mahalanobisDistSquared_x];
            obj.residuals_y = [obj.residuals_y, mahalanobisDistSquared_y];
        
            % Check for outliers using separate checks for x and y
            % Outlier rejection for x component
            if size(obj.residuals_x, 2) > M
                residual_mean_x = mean(obj.residuals_x(end-M:end-1));
                residual_std_x = std(obj.residuals_x(end-M:end-1));
        
                if mahalanobisDistSquared_x > 1 && mahalanobisDistSquared_x > (residual_mean_x + residual_std_x)
                    r_adapt_x = alpha * r_adapt_x + (1 - alpha) * (mahalanobisDistSquared_x + obj.H(1,:) * obj.P * obj.H(1,:).');
                    obj.residuals_x = obj.residuals_x(end-M:end-1);  % Keep the moving window
                end
            end
        
            % Outlier rejection for y component
            if size(obj.residuals_y, 2) > M
                residual_mean_y = mean(obj.residuals_y(end-M:end-1));
                residual_std_y = std(obj.residuals_y(end-M:end-1));
        
                if mahalanobisDistSquared_y > 1 && mahalanobisDistSquared_y > (residual_mean_y + residual_std_y)
                    r_adapt_y = alpha * r_adapt_y + (1 - alpha) * (mahalanobisDistSquared_y + obj.H(2,:) * obj.P * obj.H(2,:).');
                    obj.residuals_y = obj.residuals_y(end-M:end-1);  % Keep the moving window
                end
            end
        
            % Update the R matrix with the adapted values for x and y
            r_adapt = diag([r_adapt_x, r_adapt_y]);
        
            % Update Kalman gain with adapted covariance matrix
            S1 = obj.H * obj.P * obj.H.' + r_adapt;
            K = (obj.P * obj.H.') / S1;
        
            % Update state estimate
            obj.X = obj.X + K * (z - obj.H * obj.X);
        
            % Identity matrix for covariance update
            I = eye(size(obj.H, 2));
        
            % Update error covariance matrix
            obj.P = (I - K * obj.H) * obj.P;
            
            % Output estimated state and updated object
            X_est = obj.X;
            KF_obj2 = obj;
        end
    end

 end
        
    


