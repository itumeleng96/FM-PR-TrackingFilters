classdef GaussNewton

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,max_iter,tolerance;
    end
    
    methods
        function obj = GaussNewton(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial,max_iter,tolerance)
        
            %Init funtion
            %Inputs: 
            % dt : sampling time
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

            %Process Noise Covariance Matrix 
            obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2, 0,0,0;
                     (dt^3)/2, dt^2, dt, 0,0,0;
                     (dt^2)/2, dt, 1, 0,0,0;
                     0, 0, 0, (dt^4)/4, (dt^3)/2, (dt^2)/2;
                     0, 0, 0, (dt^3)/2, dt^2, dt;
                     0, 0, 0, (dt^2)/2, dt, 1];

            obj.Q(1:3,1:3) = obj.Q(1:3,1:3) * std_acc(1)^2;
            obj.Q(4:6,4:6) = obj.Q(4:6,4:6) * std_acc(2)^2;
            
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
            obj.X = obj.A * obj.X ;
        
            % Update the covariance matrix based on process noise
            obj.P = obj.A * obj.P * obj.A' + obj.Q;
        
            X_pred = obj.X;
            GN_Obj  = obj;

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

            lambda = 0.1;

            for i = 1:obj.max_iter
                r = obj.objectiveFunction(obj.X, z);  % Compute residual
                HtH = obj.H' * obj.H;  % Compute approximation of inverse of Jacobian
        
                while true
        
                    dx = (HtH + lambda * eye(size(HtH))) \ (obj.H' * r);
        
                    % Update state estimate and covariance matrix
                    x_new = obj.X + dx;
                    P_new = obj.A * obj.P * obj.A' + obj.Q;
                    K = P_new * obj.H' / (obj.H * P_new * obj.H' + obj.R);
                    obj.X = x_new + K * (z - obj.H * x_new);
                    obj.P = (eye(size(obj.A, 2)) - K * obj.H) * P_new;
        
                    r_new = z - obj.H * obj.X;  % Compute new residual
                    if norm(r_new) < norm(r)  % Check if iteration improved objective function
                        break;
                    end
                    lambda = lambda * 10;  % Increase damping factor
                end
                lambda = lambda * 0.1;  % Decrease damping factor
                % Check for convergence
                if norm(dx) < obj.tolerance
                    break;
                end
            end
            X_est = obj.X;
            GN_obj = obj;
        end
    end
 end
        
    


