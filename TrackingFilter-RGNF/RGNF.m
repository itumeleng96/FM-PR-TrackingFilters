classdef RGNF

    properties
        dt,U,X,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,k_d,count;
    end
    
    methods
        function obj = RGNF(dt,std_acc,r_std,rdot_std,X_initial,max_iter)
            % Constructor function to initialize the filter
            
            % Variables for Stopping criterion
            obj.max_iter = max_iter;
            
    
            % Initial State
            obj.X = X_initial;
         
            % State transition matrix
            obj.dt = dt;
            %wave number k=-lambda=c/f
            obj.k_d = -299792458/94e6; 

            %State transition matrix
            obj.A = [1,obj.k_d*dt,obj.k_d*(1/2)*dt^2;
                     0, 1, dt;
                     0, 0, 1;];
                    
            
            obj.H = [1,0,0;0,1,0;];                        % Measurement Function

            

           
            obj.Q = [(dt^4)/4,0,0;
                     0, dt^2,0;
                     0, 0, 1]*std_acc;


            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            obj.P = eye(size(obj.A,2));                    % Filter Covariance matrix
            obj.count =0;

        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
           
            obj.X = obj.A*obj.X ;
        
            % Initial covariance matrix
            obj.P = (obj.A * obj.P * obj.A.') + obj.Q;

            X_pred = obj.X;
            GN_Obj = obj;
        end


        function [X_est, RGNF_obj] = update(obj, Y_n)
            x_new = obj.X;

            %Convergence tolerance 
            tolerance = 1e-6; 
            %The forgetting  factor(Lambda) - between 0 and 1
            lambda = 1;

            %S = H*P*H'+ R - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            y = abs(x_new(2,1)-Y_n(2));


            if(y < obj.R(2,2) && obj.count<4)
                %decrease Q
                obj.count = obj.count+1;
                obj.Q = obj.Q*0.1;
            end


            for i = 1:obj.max_iter
                % Observer gain Kn
                K_n = obj.P * obj.H.' / (obj.R + obj.H * obj.P * obj.H.');
                
                % Delta Error 
                dx = K_n * (Y_n - obj.H *x_new -obj.H*(obj.X-x_new)); 
                x_temp = obj.X + dx;

                % Check for convergence
                if norm(x_temp - x_new) < tolerance
                    x_new = x_temp;
                    break;  % Converged, exit the loop
                else
                    x_new = x_temp;
                end
            end
            
            % Update Covariance Matrix
            I = eye(size(obj.H, 2));
            obj.P = (I - (K_n * obj.H)) * obj.P / lambda ;
        
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
 end