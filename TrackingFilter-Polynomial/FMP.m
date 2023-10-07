classdef FMP

    properties
        dt,U,X,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,k_d;
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

            

            obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2;
                     (dt^3)/2, dt^2, dt;
                     (dt^2)/2, dt , 1]*std_acc;


            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            obj.P = eye(size(obj.A,2));                    % Filter Covariance matrix
            obj.S = [0,0;0,0.0];                            


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
       
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
 end