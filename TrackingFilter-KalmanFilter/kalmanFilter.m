classdef kalmanFilter

    properties
        dt,U,X,F,B,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc;
    end
    
    methods
        function obj = kalmanFilter(dt,std_acc,x_std_meas,y_std_meas,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % x_std_meas : standard deviation of the measurement in the x-direction
            
            
            %Control Input Variables
            obj.U = [u_x;
                     u_y];
    
            %Initial State
            obj.X= X_initial;
                
            %Update Interval
            obj.dt = dt;

            %State transition matrix
            obj.F = [1,dt,(1/2)*dt^2;
                     0, 1, dt;
                     0, 0, 1;];


            %The control input matrix B
            obj.B = [0,0;
                     0,0;
                     0,0;];

            %Measurement Mapping Matrix 
            obj.H = [1,0,0;
                     0,1,0;];

            %Process Noise Covariance Matrix
            obj.Q = [2500, 0, 0;
                    0,1, 0;
                    0, 0,0.1];
            
            %obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2;
            %       (dt^3)/2, dt^2, dt;
            %         (dt^2)/2, dt, 1;]*100;

            %obj.Q(1:3,1:3) = obj.Q(1:3,1:3) * std_acc(1)^2;

            %Initial Measurement Noise Covariance Matrix
            %Standard deviation of measurement in doppler shift and delay
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            obj.S = [0,0;
                     0,0];
            %Initial covariance Matrix
            %High estimate uncertainty
            obj.P = eye(size(obj.F,2));

        end
        
        function [X_pred,KF_obj1] = predict(obj)
            %Prediction stage
            
            %PREDICT NEXT STATE
            %x_k = A*x_(k-1) + B*u_(k-1)
            obj.X = obj.F*obj.X ;
                        
            %COVARIANCE UPDATE
            %P= A*P*A' + Q             
            obj.P = (obj.F * obj.P) * obj.F.' + obj.Q;

            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %Update stage

            %COMPUTE KALMAN GAIN 
            
            %K = P * H'* inv(H*P*H'+R)
            %S = H*P*H'+ R - Total Error - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            K = (obj.P * obj.H.') / obj.S;
            %GET KALMAN ESTIMATE 
            obj.X = obj.X + K * (z-obj.H * obj.X);
            
            I = eye(size(obj.H,2));

            %Update Error Covariance matrix
            obj.P = (I - (K * obj.H)) * obj.P;
            
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


