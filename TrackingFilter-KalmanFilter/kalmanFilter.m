classdef kalmanFilter

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc;
    end
    
    methods
        function obj = kalmanFilter(dt,std_acc,r_std,rdot_std,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % x_std_meas : standard deviation of the measurement in the x-direction
            
        
            obj.X= X_initial;                              % Initial State
            obj.dt = dt;                                   % Update Interval
            c=299792458;    
            k = -c/94e6;                                    % Wave number k=-lambda=c/f

                    
            obj.F = [1, 0, k*dt,0;
                     0, 0, k, k*dt;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];                             % Measurement Function

            

            obj.Q = [5,0,0,0;
                     0, 0.02, 0, 0;
                     0, 0, 0.2,0;
                     0, 0, 0, 0.05];


            obj.R = [r_std,0;0,rdot_std];                  % Measurement Uncertainty
            
            obj.P = [5,0,0,0;                              % Initial Error Covariance Matrix
                     0, 0.02, 0, 0;
                     0, 0, 0.4,0;
                     0, 0, 0, 0.1];                           

            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];


        end
        
        function [X_pred,KF_obj1] = predict(obj)

            %Predict next state (prior)
           
            % x = Fx
            obj.X = obj.F*obj.X ; 
            
            % P = FPF' + Q
            obj.P = obj.A * obj.P * obj.A.' + obj.Q;

            %Return Prior
            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %UPDATE STAGE
            
            % Doppler measurement is reliable, perform the Kalman update
            %S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;


            %K = PH'inv(S)
            
            K = (obj.P * obj.H.') * obj.S^(-1);
            %x = x + R
            obj.X = obj.X + K * (z-obj.H * obj.X);
            I = eye(size(obj.H,2));
    
            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P ;            
          
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


