classdef EKF

    properties
        dt,U,X,A,B,H,Q,R,P;
    end
    
    methods
        function obj = EKF(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas)
        
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
            obj.X= [16;
                    0;
                    370;
                    0];

            %State transition matrix
            obj.A = [1,dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1];
                 
             

            %The control input matrix B
            obj.B = [(dt^2)/2, 0;
                     0, (dt^2)/2;
                     (dt^2)/2, 0;
                     0, (dt^2)/2;];

            %Measurement Mapping Matrix -must compute Jacobian
            obj.H = [1,0,0,0;
                     0,0,1,0];

            %Process Noise Covariance
            %obj.Q = [(dt^4)/4  , 0, (dt^3)/2, 0;
            %         0, (dt^4)/4, 0, (dt^3)/2;
            %         (dt^3)/2, 0, dt^2, 0;
            %         0, (dt^3)/2, 0, dt^2] .*std_acc^2;
            obj.Q = [1,0,0,0;
                     0,1,0,0;
                     0,0,1,0;
                     0,0,0,1].*std_acc^2;
            
            %obj.Q = [(dt^4)/2, 0, dt^3, 0;
            %         0, (dt^4)/2, 0, dt^3;
            %        (dt^4)/2, 0, dt^3, 0;
            %         0, (dt^4)/2, 0, dt^3] .*std_acc^2;

            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                 0, y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.A,2));
        end
        
        function [X_pred,KF_obj1] = predict(obj)
            %Calculate the predicted time state

            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            obj.X= obj.A * obj.X + obj.B * obj.U;
            %calculate error covariance
            %P= A*P*A' + Q 
            obj.P = (obj.A * obj.P) * obj.A.' + obj.Q;
            
            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %Update stage, compute the Kalman gain

            %S = H*P*H'+R
            S = obj.H * obj.P * obj.H.' + obj.R ;

            %Calculate the Kalman Gain 
            %K = P * H'* inv(H*P*H'+R)

            K = (obj.P * obj.H.') * S^(-1);
            
            obj.X = round(obj.X + K * (z-obj.H * obj.X));
                       
            I = eye(size(obj.H,2));

            %Update Error Covariance matrix
            obj.P = (I - (K * obj.H)) * obj.P;
            
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
        function [Aj,Hj]  = recompute(obj)
            Aj = jacobian(obj.A,obj.X);
            Hj = jacobian(obj.H,obj.X);
        end

     end
 end
        
    


