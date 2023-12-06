classdef FMP

    properties
        dt,U,X,A,B,H,R,S,measured_x,measured_y,n,k_d;
    end
    
    methods
        function obj = FMP(dt,std_acc,r_std,rdot_std,X_initial)
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
           

            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            obj.S = [0,0;0,0.0];  
            obj.n =0;

        end

        % Function to predict the next state
        function [X_pred, FMP_Obj] = predict(obj)
           
            obj.X = obj.A*obj.X ;
        
            X_pred = obj.X;
            FMP_Obj = obj;
        end


        function [X_est, FMP_obj] = update(obj, Y_n)
            %2nd degree Fading Memory Polynomial Filter
            
            X_k = obj.X;
            %Calculate the Error
            e_n = Y_n-obj.H*X_k;
            
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

            Temp0=z0+alpha*e_n(1);
            Temp1=z1+beta*e_n(1);
            
            Temp2=z2+alpha*e_n(2);
            Temp3=z3+beta*e_n(2);
            
            
            
            Z=[Temp0;Temp1;Temp2;Temp3];

            obj.X = Z;
            %Update Batch Number or Update number
            obj.n = obj.n+1;
            X_est = obj.X;
            FMP_obj = obj;

        end

    end
 end