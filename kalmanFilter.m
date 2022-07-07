classdef kalmanFilter 

    properties
        dt,U,X,A,B,H,Q,R,P;
    end
    
    methods
        function obj = kalmanFilter(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % u_x : acceleration in the x-direction
            % u_y : acceleration in the y-direction
            % x_std_meas : standard deviation of the measurement in the x-direction
            % y_std_meas : standard deviation of the measurement in the y-direction
            
            %Control Input Variables
            obj.U = [u_x,u_y];

            %Initial State
            obj.X= zeros(1,4);

            %State transition matrix
            obj.A = [1, 0,dt, 0;
                 0, 1, 0,dt;
                 0, 0, 1, 0;
                 0, 0, 0, 1];

            %The control input matrix B
            obj.B = [(dt^2)/2, 0;
                 0, (dt^2)/2;
                 dt, 0;
                 0,dt];

            %Measurement Mapping Matrix 
            obj.H = [1,0,0,0;
                 0,1,0,0];

            %Process Noise Covariance
            obj.Q = [(dt^4)/4  , 0, (dt^3)/2, 0;
                 0, (dt^4)/4, 0, (dt^3)/2;
                 (dt^3)/2, 0, dt^2, 0;
                 0, (dt^3)/2, 0, dt^2] .*std_acc^2;

            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                 0, y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(A));
        end
        
        function X = predict(obj)
            %Calculate the predicted time state

            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            obj.X= dot(obj.A,obj.X) + dot(obj.B,obj.U);

            %calculate error covariance
            %P= A*P*A' + Q 
            obj.P = dot(dot(obj.A,obj.P),obj.A.') + obj.Q;
            
            X = obj.X;
        end
        
        function X = update(obj,z)
            %Update stage, compute the Kalman gain

            %S = H*P*H'+R
            S = dot(obj.H,dot(obj.P,obj.H.')) + obj.R ;

            %Calculate the Kalman Gain 
            %K = P * H'* inv(H*P*H'+R)

            K = dot(dot(obj.P,obj.H.'),S^(-1));

            obj.X = round(obj.X+ dot(K,(z-dot(obj.H,obj.X))));

            I = eye(size(obj.H));

            %Update Error Covariance matrix
            obj.P = (I - (K*obj.H))* obj.P;
            
            X = obj.X;
        end

     end
 end
        
    


