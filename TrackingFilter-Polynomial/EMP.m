classdef EMP

    properties
        dt,U,X,A,B,H,R,S,measured_x,measured_y,n,k_d;
    end
    
    methods
        function obj = EMP(dt,std_acc,r_std,rdot_std,X_initial)
            % Constructor function to initialize the filter
            
    
            % Initial State
            obj.X = X_initial;
         
            % State transition matrix
            obj.dt = dt;
            %wave number k=-lambda=c/f
            obj.k_d = -299792458/94e6; 

            %State transition matrix
            obj.A = [1,dt,(1/2)*dt^2;
                     0, 1, dt;
                     0, 0, 1;];
                    
            
            obj.H = [1,0,0;0,1,0;];                        % Measurement Function

           

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
            %Calculate the Error
            e_n = [Y_n;0] -X_k;
            disp(e_n);
            %Calculate weights
            n = obj.n;
            gamma = 30 / ((n + 3)*(n+2)*(n+1));
            beta = 18 * (2 * n + 1) / ((n + 3)*(n+2)*(n+1));
            alpha = 3 * (3 * n^2 + 3 * n + 2) /((n + 3)*(n+2)*(n+1));
            %T_n = [alpha;beta;gamma;];

            %Update
            z0=obj.X(1);
            z1=obj.X(2);
            z2=obj.X(3);

            Temp0=z0+z1+z2+alpha*e_n(1);
            Temp1=z1+2*z2+beta*e_n(2);
            Temp2=z2+gamma*e_n(2);
            
            
            Z=[Temp0;Temp1;Temp2];

            disp(Z);
            obj.X = Z;
            %Update Batch Number or Update number
            obj.n = obj.n+1;
            X_est = obj.X;
            EMP_obj = obj;

        end

    end
 end