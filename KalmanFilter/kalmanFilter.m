classdef kalmanFilter

    properties
        dt,U,X,F,B,H,Q,R,P,coeff,measured_x,measured_y;
    end
    
    methods
        function obj = kalmanFilter(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % u_x : acceleration in the x-direction
            % u_y : acceleration in the y-direction
            % x_std_meas : standard deviation of the measurement in the x-direction
            % y_std_meas : standard deviation of the measurement in the y-direction
            
            
            %Control Input Variables
            obj.U = [u_x;
                     u_y];
    
            %Initial State
            obj.X= X_initial;
         

            %State transition matrix
            obj.dt = dt;
            lambda = 1e-9;

            obj.F = [1,0,-lambda*dt,0;
                     0, 0,-lambda,-lambda*dt;
                     0, 0, 1, dt;
                     0, 0  0, 1;];
                    
            obj.H = [1,0,0,0;
                     0,0,1,0];

            obj.Q = [0.01,0,0,0;
                     0, 0.02, 0,0;
                     0, 0, 0.2, 0;
                     0, 0  0, 0.05;];

            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.F,2));

        end
        
        function [X_pred,KF_obj1] = predict(obj)
            %Calculate the predicted time state
            
            %Update time state
            %x_k = F.x_(k-1) + B.u_(k-1) 
            obj.X= obj.F * obj.X ;
                        
            %calculate error covariance
            %P= F*P*F' + Q 
            obj.P = eye(size(obj.F,2));
            
            obj.P = (obj.F * obj.P) * obj.F.' + obj.Q;
            
            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %Update stage, compute the Kalman gain

            %S = H*P*H'+ R - Total Error 
            S = obj.H * obj.P * obj.H.' + obj.R;

            %Calculate the Kalman Gain 
            %K = P * H'* inv(H*P*H'+R)

            K = (obj.P * obj.H.') * S^(-1);
            
            obj.X = obj.X + K * (z-obj.H * obj.X);
                        
            I = eye(size(obj.H,2));

            %Update Error Covariance matrix
            obj.P = (I - (K * obj.H)) * obj.P;
            
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


