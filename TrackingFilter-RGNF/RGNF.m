classdef RGNF

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,max_iter,tolerance;
    end
    
    methods
        function obj = RGNF(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial,max_iter,tolerance)
            % Constructor function to initialize the filter
            
            % Variables for Stopping criterion
            obj.max_iter = max_iter;
            obj.tolerance = tolerance;
            
            % Control Input Variables
            obj.U = [u_x; u_y];
    
            % Initial State
            obj.X = X_initial;
         
            % State transition matrix
            obj.dt = dt;
            obj.A = [1, dt, (1/2)*dt^2;
                     0, 1, dt;
                     0, 0, 1;];
                 
            % The control input matrix B
            obj.B = [0, 0;
                     0, 0;
                     0, 0;];
                 
            % Measurement Mapping Matrix 
             obj.H = [1,0,0;
                     0,1,0;];
                 
            obj.Q = [100, 0, 0;
                    0,1, 0;
                    0, 0,0.1];
            
            %obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2;
            %       (dt^3)/2, dt^2, dt;
            %         (dt^2)/2, dt, 1;]*100;
                 
            %obj.Q(1:3, 1:3) = obj.Q(1:3, 1:3) * std_acc(1)^2;
            
            % Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2, 0;
                     0, y_std_meas^2];
        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
            % Calculate the predicted time state: x_k = Ax_(k-1) + Bu_(k-1)
            % Generate random Gaussian noise with zero mean and covariance matrix Q
        
            % Generate random Gaussian noise with zero mean and covariance matrix Q
            noise = sqrtm(obj.Q) * randn(size(obj.X));


            obj.X = obj.A * obj.X +noise;
        
            % Update the covariance matrix based on process noise
            X_pred = obj.X;
            GN_Obj = obj;
        end

        function [residual] = objectiveFunction(obj, X_pred, z)
            % Compute the residual: residual = z - (H * X_pred)
            residual = z - (obj.H * X_pred);
        end

        function [X_est, RGNF_obj] = update(obj, z)
            % Perform the update step of the Gauss-Newton filter

            lambda = 0.1;

            for i = 1:obj.max_iter
                r = obj.objectiveFunction(obj.X, z);  % Compute residual
                HtH = obj.H' * obj.H;                 % Compute approximation of inverse of Jacobian
        
                while true
                    dx = (HtH + lambda * eye(size(HtH))) \ (obj.H' * r);
        
                    % Update state estimate and covariance matrix
                    x_new = obj.X + dx;
                    obj.X = x_new;        
                    r_new = z - obj.H * obj.X;  % Compute new residual
                    if norm(r_new) < norm(r)  % Check if iteration improved objective function
                        break;
                    end
                    lambda = lambda * 10;  % Increase damping factor
                end
                
                lambda = lambda * 0.1;     % Decrease damping factor
                
                % Check for convergence
                if norm(dx) < obj.tolerance
                    break;
                end
            end
            
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
 end