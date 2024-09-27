classdef RGNF

    properties
        dt,U,X,F,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,wk,count,updater,update1,lambda;
    end
    
    methods
        function obj = RGNF(dt,std_acc,r_std,rdot_std,X_initial,max_iter,lambda)
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

            

             obj.Q = [std_acc(1)*(dt^4)/4,std_acc(1)*(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                 std_acc(1)*(dt^3)/2, std_acc(1)*dt^2, 0, 0;
                 0, 0, std_acc(2)*(dt^4)/4,std_acc(2)*(dt^3)/2;
                 0, 0, std_acc(2)*(dt^3)/2, std_acc(2)*dt^2];


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
            
            obj.wk = [std_acc(1)*dt^2;std_acc(1)*dt;std_acc(2)*dt^2;std_acc(2)*dt];
            obj.lambda = lambda;

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
            obj.count = obj.count+1;    

            %The forgetting  factor(Lambda) - between 0 and 1

            %S = H*P*H'+ R - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;


            for i = 1:obj.max_iter
                % Observer gain Kn
                K_n = obj.P * obj.H.' *obj.S^(-1);
                
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

            obj.P = (I - (K_n * obj.H)) * obj.P / obj.lambda ;
        
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
 end