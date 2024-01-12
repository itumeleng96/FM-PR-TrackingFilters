classdef EMP

    properties
        dt,dp,U,X,A,B,H,R,S,measured_x,measured_y,n,k_d;
    end
    
    methods
        function obj = EMP(dt,std_acc,r_std,rdot_std,X_initial)
            % Constructor function to initialize the filter
            
    
            % Initial State
            obj.X = X_initial;
         
            % State transition matrix
            obj.dt = dt;

            %State transition matrix
            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];     

            obj.dp = [1;dt;1;dt;];
            

            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            obj.S = [0,0;0,0.0];  
            obj.n =0;

        end

        % Function to predict the next state
        function [X_pred, EMP_Obj] = predict(obj)
           
            obj.X = obj.A*obj.X ;
        

            X_pred = obj.X;
            EMP_Obj = obj;
        end


        function [X_est, EMP_obj] = update(obj, Y_n)
            %2nd degree Expanding Memory Polynomial Filter
            
            X_k = obj.X;
            Z_k = X_k.*obj.dp;

            %Calculate the Error
            E_n = Y_n-obj.H*X_k;
            
            %Calculate weights
            beta = 6/ ((obj.n + 2)*(obj.n+1));
            alpha = 2 * (2 * obj.n + 1) /((obj.n + 2)*(obj.n+1));
            

            Temp3= Z_k(4) + beta *E_n(2);
            Temp2 = Z_k(3) + alpha* E_n(2);
            Temp1 = Z_k(2) + beta *E_n(1);
            Temp0 = Z_k(1) +  alpha * E_n(1);
   
            
            X_k = [Temp0;(1/obj.dt)*Temp1;Temp2;(1/obj.dt)*Temp3];

            obj.X = X_k;
            X_est = obj.X;

            obj.n = obj.n+1;
            
            EMP_obj = obj;
        end

    end
 end