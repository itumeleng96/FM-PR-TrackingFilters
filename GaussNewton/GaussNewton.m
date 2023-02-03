classdef GaussNewton

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,max_iter,tolerance;
    end
    
    methods
        function obj = GaussNewton(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial,max_iter,tolerance)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % u_x : acceleration in the x-direction
            % u_y : acceleration in the y-direction
            % x_std_meas : standard deviation of the measurement in the x-direction
            % y_std_meas : standard deviation of the measurement in the y-direction
            
            %variables for Stopping criterion
            obj.max_iter=max_iter;
            obj.tolerance = tolerance;
            
            %Control Input Variables
            obj.U = [u_x;
                     u_y];
    
            %Initial State
            obj.X= X_initial;
         

            %State transition matrix
            obj.dt = dt;


            obj.A = [1,dt,(1/2)*dt^2, 0 ,0 ,0;
                     0, 1, dt,0,0,0;
                     0, 0, 1, 0,0,0;
                     0, 0  0, 1,dt,(1/2)*dt^2;
                     0, 0, 0, 0,1,dt ;
                     0, 0, 0, 0,0,1;];


            %The control input matrix B
            obj.B = [0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;];

            %Measurement Mapping Matrix 
            obj.H = [1,0,0,0,0,0;
                     0,0,0,1,0,0];

            %Process Noise Covariance
            obj.Q = [0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0];

            
            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.A,2));

        end

        %Function to predict the next state
        function [X_pred,GN_Obj] = predict(obj)
            %Calculate the predicted time state
            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            obj.X= obj.A * obj.X + obj.B * obj.U;
                        
            
            X_pred = obj.X;
            GN_Obj  = obj;
            disp("Prediction");
            disp(X_pred);
        end

        %sum of squared differences
        function [residual] = objectiveFunction(obj,X_pred,z)
            %residual = (z-(obj.H*X_pred))' * (z-(obj.H*X_pred));
            residual= z-(obj.H*X_pred);
        end
        %{
        function [J] =jacobian(obj)
            %J_h(i, j) = ∂h(p)/∂p_j = ∂h_i/∂p_j

            % J = -H./(H*x + R);
            % Jacobian matrix
            J = -obj.H./(obj.H*obj.X);

        end
        %}

        function [X_est,GN_obj] = update(obj,z)
            lambda = 0.01;

            for i = 1:obj.max_iter
                r = obj.objectiveFunction(obj.X,z);         % Compute residual
                HtH = obj.H'*obj.H;                         % Compute approximation of inverse of Jacobian

                while true

                    dx = (HtH + lambda * eye(size(HtH))) \ (obj.H' * r);
                    
                    x_new = obj.X + dx;                                % Update state estimate
                    r_new = z - obj.H*x_new;                           % Compute new residual
                    if norm(r_new) < norm(r)                           % Check if iteration improved objective function
                        obj.X = x_new;                                     % Accept new state estimate
                        break;
                    end
                    lambda = lambda * 10; % Increase damping factor
                end
                lambda = lambda * 0.1; % Decrease damping factor
                % Check for convergence
                if norm(dx) < 1e-6
                    break;
                end
            end 
            X_est = obj.X;
            GN_obj=obj;
        end
     end
 end
        
    


