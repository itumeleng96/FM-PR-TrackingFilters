classdef kalmanFilter

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,k_d,w_k,count;
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
            obj.k_d = -c/94e6;                             % Wave number k=-lambda=c/f

                    
            obj.F = [1,obj.k_d*dt,dt^2;
                     0, 1, dt;
                     0, 0, 1;];
                    
            
            obj.H = [1,0,0;0,1,0;];                        % Measurement Function

            

            obj.Q = [1500,0,0;
                     0, 1,0;
                     0, 0, 0.1]*std_acc;


            obj.R = [r_std^2,0;0,rdot_std^2];              % Measurement Uncertainty
            
            obj.P = eye(size(obj.F,2));                    % Uncertainty Covariance
           
            obj.A = [1, dt, dt^2;
                     0, 1, dt ;
                     0, 0, 1];

            obj.w_k = [(dt^2)/2;dt;1];
            obj.count =0;
      
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
            
            %y = abs(obj.X(2,1)-z(2));

            %if(y < obj.R(2,2) && obj.count<4)
                %decrease Q
             %   obj.count = obj.count+1;
             %   obj.Q = obj.Q*0.1;
            %end


            %K = PH'inv(S)
            K = (obj.P * obj.H.') * obj.S^(-1);
            %x = x + Ky
            obj.X = obj.X + K * (z-obj.H * obj.X);1
            I = eye(size(obj.H,2));
    
            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P ;

            
          
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


