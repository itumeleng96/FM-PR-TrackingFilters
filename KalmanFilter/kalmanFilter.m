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


            obj.F = [1,dt,(1/2)*dt^2, 0 ,0 ,0;
                     0, 1, dt,0,0,0;
                     0, 0, 1, 0,0,0;
                     0, 0  0, 1,dt,(1/2)*dt^2;
                     0, 0, 0, 0,1,dt ;
                     0, 0, 0, 0,0,1;];


            %The control input matrix B
            obj.B = [0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;];

            %Measurement Mapping Matrix 
            obj.H = [1,0,0,0,0,0;
                     0,0,0,1,0,0];

            %Process Noise Covariance
            obj.Q = [(dt^5)/20, (dt^4)/8, (dt^3)/6, 0,0,0;
                     (dt^4)/8, dt^3/3, dt^2/2, 0,0,0;
                     (dt^3)/6, dt^2/2, dt, 0,0,0;
                     0, 0, 0, (dt^5)/20, (dt^4)/8, (dt^3)/6;
                     0, 0, 0, (dt^4)/8, dt^3/3, dt^2/2;
                     0, 0, 0, (dt^3)/6, dt^2/2, dt].*0.5;

            
            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.F,2));

        end
        
        function [X_pred,KF_obj1] = predict(obj)
            %Calculate the predicted time state
            
            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            obj.X= obj.F * obj.X ;
                        
            %calculate error covariance
            %P= A*P*A' + Q 
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
        
    


