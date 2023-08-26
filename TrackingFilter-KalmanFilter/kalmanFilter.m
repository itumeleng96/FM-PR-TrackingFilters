classdef kalmanFilter

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,k_d;
    end
    
    methods
        function obj = kalmanFilter(dt,std_acc,r_std,rdot_std,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % x_std_meas : standard deviation of the measurement in the x-direction
            
        
            %Initial State
            obj.X= X_initial;
                
            %Update Interval
            obj.dt = dt;

            %wave number k=-lambda=c/f
            obj.k_d = -299792458/94e6; 

            %State transition matrix
            obj.F = [1, obj.k_d*dt;
                     0, 1];

            %Transition Matrix for P(Error-Covariance)
            obj.A = [1, dt;
                     0, 1];

            %Measurement Mapping Matrix 
            obj.H = [1,0;
                     0,1;];

            %Process Noise Covariance Matrix For Random Acceleration
            obj.Q = [(dt^4)/4,(dt^3)/2;
                    (dt^3)/2,dt^2]*std_acc;
            
            %Standard deviation of measurement in doppler shift and delay
            %Measurement Error covariance matrix
            obj.R = [r_std^2,0.01*rdot_std*r_std;
                     0.01*rdot_std*r_std,rdot_std^2];

            %Initial Innovation Error Matrix
            obj.S = [0,0;
                     0,0];

            obj.P = eye(size(obj.A,2));

        end
        
        function [X_pred,KF_obj1] = predict(obj)
            %PREDICTION STAGE
            
            %x_k = F*x_(k-1)
            obj.X = obj.F*obj.X ;
                       
            %P= A*P*A' + Q             
            obj.P = obj.A * obj.P * obj.A.' + obj.Q;

            X_pred = obj.X;
            KF_obj1  = obj;
        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %UPDATE STAGE
            
            %S = H*P*H'+ R - Total Error - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            %KALMAN GAIN
            %K = P * H'* inv(H*P*H'+R)
            K = obj.P * obj.H.'*obj.S^(-1);

            %GET ESTIMATE 
            obj.X = obj.X + K * (z-obj.H * obj.X);
            
            I = eye(size(obj.H,2));

            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P;
            
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


