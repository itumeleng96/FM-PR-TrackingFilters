classdef CSRGNF

    properties
        dt,U,X,F,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,wk,count,updater,update1,residuals;
    end
    
    methods
        function obj = CSRGNF(dt,std_acc,r_std,rdot_std,X_initial,max_iter)
            % Constructor function to initialize the filter
            
            % Variables for Stopping criterion
            obj.max_iter = max_iter;
            
    
            % Initial State
            obj.X = X_initial;
         
            obj.dt = dt;

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

            %Measurement Error covariance matrix
            obj.R = [r_std,0;
                     0,rdot_std;];

            obj.P = [5,0,0,0;                                      % Initial Error Covariance Matrix
                     0, 2.5, 0, 0;
                     0, 0, 2,0;
                     0, 0, 0, 1];  


            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
            
            obj.wk = std_acc*[dt^2;dt;dt^2;dt];

            obj.residuals =[];

        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
           
            obj.X = obj.F*obj.X;
        
            % Initial covariance matrix
            obj.P = (obj.A * obj.P * obj.A.') + obj.Q;
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            X_pred = obj.X;
            GN_Obj = obj;
        end


        function [X_est, RGNF_obj] = update(obj, Y_n)
            x_new = obj.X;

            %Convergence tolerance 
            tolerance = 1e-1; 

            %The forgetting  factor(Lambda) - between 0 and 1
            lambda = 1;

            % Kalman prediction step: S = H*P*H' + R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            
            % Compute residual (measurement innovation): ek = z - H*X
            ek = Y_n - obj.H * obj.X;
            
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Recursive Gauss Newton Update Step
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           

            S1 = obj.H * obj.P * obj.H.' + r_adapt;
            I = eye(size(obj.H, 2));

            % Recursive Update
            for i = 1:obj.max_iter                
                % Observer gain Kn
                % K_n = Wn^-1_n-1*H^T[R_n + H W^-1_n-1 H^T]^-1
                K_n = obj.P * obj.H.' *S1^(-1); 
                dx = K_n*(Y_n - obj.H * x_new - obj.H * (obj.X - x_new));
                %weights = obj.huberFunction(obj.R,residual, 1.5); % Compute weights using Huber loss function
            
                % Compute the update
                x_temp = obj.X + dx;

                % Check for convergence
                if norm(x_temp-x_new) < tolerance
                    x_new = x_temp;
                    % Converged, exit the loop
                    break; 
                else
                    x_new = x_temp;
                end
            end

            % W^−1_n = λ^−1[I −Kn H]W^−1_n-1
            obj.P = (I - K_n * obj.H) * obj.P / lambda;

                            
            % Update Covariance Matrix
        
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
    methods(Static)
        function weights = huberFunction(residual, delta)                    
            weights = zeros(size(residual));
            quadratic_indices = abs(norm(res)) <= delta;
            
            weights(quadratic_indices) = 1; % Quadratic loss
            weights(~quadratic_indices) = delta ./ abs(residual(~quadratic_indices)); % Linear loss
        end
    end
    
 end
