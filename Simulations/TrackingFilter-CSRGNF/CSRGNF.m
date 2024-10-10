classdef CSRGNF

    properties
        dt,U,X,F,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,wk,count,updater,update1,epsDoppler,epsRange;
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

            obj.epsDoppler =[];
            obj.epsRange =[];

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
        
            % Convergence tolerance 
            tolerance = 1e-1; 
        
            % The forgetting factor (Lambda) - between 0 and 1
            lambda = 1;
        
            % Kalman prediction step: S = H*P*H' + R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
        
            % Compute residual (measurement innovation) for both dimensions
            ek = Y_n - obj.H * obj.X;
            
            % Mahalanobis distances for x (Range) and y (Doppler)
            eps_range = (ek(1)^2) / obj.S(1,1);  % Range residual (x component)
            eps_doppler = (ek(2)^2) / obj.S(2,2);  % Doppler residual (y component)
        
            % Store residuals in moving window arrays
            obj.epsRange = [obj.epsRange, eps_range];
            obj.epsDoppler = [obj.epsDoppler, eps_doppler];

        
            % Number of samples to average
            M = 6;
            alphaFactor = 0.6;
        
            % Initialize adaptive noise covariance
            r_adapt = obj.R; 
        
        
            % Outlier rejection for Doppler
            if (size(obj.epsDoppler, 2) > M)
                residual_mean_doppler = mean(obj.epsDoppler(end-M:end-1));
                residual_std_doppler = std(obj.epsDoppler(end-M:end-1));
                
                if eps_doppler > 1 && eps_doppler > (residual_mean_doppler + residual_std_doppler)
                    r_adapt(2, 2) = alphaFactor * r_adapt(2, 2) + (1 - alphaFactor) * (eps_doppler + obj.P(3, 3));
                    % Remove old values from moving window
                    obj.epsDoppler = obj.epsDoppler(end-M:end-1);
                end
            end
        
            % Outlier rejection for Range (x component)
            if (size(obj.epsRange, 2) > M)
                residual_mean_range = mean(obj.epsRange(end-M:end-1));
                residual_std_range = std(obj.epsRange(end-M:end-1));
                
                if eps_range > 2 && eps_range > (residual_mean_range + residual_std_range)
                    r_adapt(1, 1) = alphaFactor * r_adapt(1, 1) + (1 - alphaFactor) * (eps_range + obj.P(1, 1));
                    % Remove old values from moving window
                    obj.epsRange = obj.epsRange(end-M:end-1);
                end
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Recursive Gauss Newton Update Step
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % Update S1 using the adapted noise covariance
            S1 = obj.H * obj.P * obj.H.' + r_adapt; 
            I = eye(size(obj.H, 2));
        
            % Recursive Update
            for i = 1:obj.max_iter                
                % Observer gain Kn
                K_n = obj.P * obj.H.' * S1^(-1); 
                dx = K_n * (Y_n - obj.H * x_new - obj.H * (obj.X - x_new));
                
                % Compute the update
                x_temp = obj.X + dx;
        
                % Check for convergence
                if norm(x_temp - x_new) < tolerance
                    x_new = x_temp;
                    % Converged, exit the loop
                    break; 
                else
                    x_new = x_temp;
                end
            end
        
            % Update the covariance matrix
            obj.P = (I - K_n * obj.H) * obj.P / lambda;
        
            % Update state estimate
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
