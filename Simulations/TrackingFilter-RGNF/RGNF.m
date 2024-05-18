classdef RGNF

    properties
        dt,U,X,F,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,wk,count,updater,update1;
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

            

            obj.Q = [(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            
            obj.P = [1,0,0,0;                              % Initial Error Covariance Matrix
                     0, 0.02, 0, 0;
                     0, 0, 0.4,0;
                     0, 0, 0, 0.1];                           

            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
            
            obj.wk = std_acc*[dt^2;dt;dt^2;dt];


            obj.count =0;
            obj.updater =0;
            obj.update1 =0;
        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
           
            obj.X = obj.F*obj.X +obj.wk ;
        
            % Initial covariance matrix
            obj.P = (obj.A * obj.P * obj.A.') + obj.Q;

            X_pred = obj.X;
            GN_Obj = obj;
        end


        function [X_est, RGNF_obj] = update(obj, Y_n)
            x_new = obj.X;

            %Convergence tolerance 
            tolerance = 1e-6; 
            threshold=3;
            max_adapt=4;
            obj.count = obj.count+1;    

            %The forgetting  factor(Lambda) - between 0 and 1
            lambda = 1;

            %S = H*P*H'+ R - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            eps_ =log(normpdf(Y_n(2),obj.X(3,1),obj.S(2,2)));
            Kp = 0;
            
            if(abs(eps_)>threshold && obj.count>10 && obj.updater<max_adapt && obj.update1>max_adapt)
                Kp = (obj.P * obj.H.') * obj.S^(-1);

                obj.updater = obj.updater+1;
                if(obj.updater==max_adapt)
                    obj.update1=0;
                end

                r_adapt = obj.R;
                r_adapt(2,2) =r_adapt(2,2)*1000;
                obj.S = obj.H * obj.P * obj.H.' + r_adapt;
            else
                obj.updater =0;
                obj.update1=obj.update1+1;
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
            if(Kp==0)
                obj.P = (I - (K_n * obj.H)) * obj.P / lambda ;

            else
                obj.P = (I - (Kp * obj.H)) * obj.P / lambda ;

            end
        
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
 end