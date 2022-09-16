classdef EKF

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y;
    end
    
    methods
        function obj = EKF(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial,coeff)
        
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

            obj.A = [1,dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1];
                 
            obj.dt = dt;
            obj.coeff  = coeff;
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
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
        end
        
        function [X_pred,KF_obj1] = predict(obj)
            %Calculate the predicted time state
            
            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            %obj.X= obj.A * obj.X + obj.B * obj.U;
            
            %For the predict function a non-linear quadratic function is
            %used
            X_ = zeros(4);
            X_(1) = obj.X(1) + obj.X(2)*obj.dt;
            X_(2) = obj.X(2);
            X_(3) = obj.coeff(1)*obj.X(1)^2 +obj.coeff(2)*obj.X(1) + obj.coeff(3);
            X_(4) = obj.X(4);
            disp(obj.coeff);

            %Calculate the Jacobian Matrix
            a31 = 2*obj.coeff(1)*obj.X(1);
            a32 = obj.coeff(2);
            A_k = [1,obj.dt, 0, 0;
                   0, 1, 0, 0;
                   a31, a32, 0,0;
                   0, 0, 0, 1];
            %calculate error covariance
            %P= A*P*A' + Q 
            obj.P = eye(size(obj.A,2));
            obj.P = (A_k * obj.P) * A_k.' + obj.Q;
            
            obj.X = X_;
            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z,x,y,index)
            %Update stage, compute the Kalman gain
            obj.measured_y =y;
            obj.measured_x =x;
            if index>1
                obj.coeff = polyfit(obj.measured_x,obj.measured_y,2);
            end

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
        
     end
 end
        
    


