classdef RGNF

    properties
        dt,U,X,F,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter;
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
            k= -299792458/94e6; 

            obj.F = [1, 0, k*dt,0;
                     0, 0, k, k*dt;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];                             % Measurement Function

            

            obj.Q = [1,0,0,0;
                     0, 0.02, 0, 0;
                     0, 0, 0.2,0;
                     0, 0, 0, 0.005];


            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            
            obj.P = [1,0,0,0;                              % Initial Error Covariance Matrix
                     0, 0.02, 0, 0;
                     0, 0, 0.4,0;
                     0, 0, 0, 0.1];                           

            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];

        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
           
            obj.X = obj.F*obj.X ;
        
            % Initial covariance matrix
            obj.P = (obj.A * obj.P * obj.A.') + obj.Q;

            X_pred = obj.X;
            GN_Obj = obj;
        end


        function [X_est, RGNF_obj] = update(obj, Y_n)
            x_new = obj.X;

            %Convergence tolerance 
            tolerance = 1e-6; 
            threshold=10;

            %The forgetting  factor(Lambda) - between 0 and 1
            lambda = 1;

            %S = H*P*H'+ R - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            dk = obj.X(3,1)-Y_n(2);
            dk_average = dk*dk;
            alpha = dk_average/obj.S(2,2);

            if(abs(alpha)>threshold)
                r_adapt = obj.R;
                r_adapt(2,2) =r_adapt(2,2)*1000;
                obj.S = obj.H * obj.P * obj.H.' + r_adapt;
            else
                obj.S = obj.H * obj.P * obj.H.' + obj.R;
            end
            


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
            obj.P = (I - (K_n * obj.H)) * obj.P / lambda ;
        
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
 end