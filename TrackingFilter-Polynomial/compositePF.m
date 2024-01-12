classdef compositePF

    properties
        dt,dp,U,X,A,B,H,R,S,measured_x,measured_y,n,k_d;
    end
    
    methods
        function obj = compositePF(dt,std_acc,r_std,rdot_std,X_initial)
            % Constructor function to initialize the filter
            
    
            % Initial State
            obj.X = X_initial;
         
            % State transition matrix
            obj.dt = dt;
            %wave number k=-lambda=c/f
            obj.k_d = -299792458/94e6; 

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
            %Composite filter
            %Initializes with the EMP filter and transitions to FMP
            
            Ns=24;

            X_k = obj.X;
            Z_k = X_k.*obj.dp;

            %Calculate the Error
            E_n = Y_n-obj.H*X_k;
            
            if (obj.n<Ns)
                %Calculate weights using EMP
                %Calculate weights
                beta = 6/ ((obj.n + 2)*(obj.n+1));
                alpha = 2 * (2 * obj.n + 1) /((obj.n + 2)*(obj.n+1));
                
                Temp3= Z_k(4) + beta *E_n(2);
                Temp2 = Z_k(3) + alpha* E_n(2);
                Temp1 = Z_k(2) + beta *E_n(1);
                Temp0 = Z_k(1) +  alpha * E_n(1);
   
            
            X_k = [Temp0;(1/obj.dt)*Temp1;Temp2;(1/obj.dt)*Temp3];

            else
                %Calculate weights
                theta = 0.83;
                beta =  (1-theta)^2;
                alpha = 1 - theta^2;
    
                %T_n = [alpha;beta;gamma;];
    
                %Update
                z0=obj.X(1);
                z1=obj.X(2);
                z2=obj.X(3);
                z3=obj.X(4);
    
                Temp0=z0+alpha*E_n(1);
                Temp1=z1+beta*E_n(1);
                
                Temp2=z2+alpha*E_n(2);
                Temp3=z3+beta*E_n(2);
                
                
                
                X_k=[Temp0;Temp1;Temp2;Temp3];

            end

            obj.X = X_k;
            X_est = obj.X;

            obj.n = obj.n+1;
            
            EMP_obj = obj;
        end

    end
 end