classdef kalmanFilter

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,dk_history,P_adapt,wk,count,counter,updater,update1;
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

            

            obj.Q = [(dt^4)/4,(dt^3)/2,0,0;                              % Initial Error Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            obj.R = [r_std,0;0,rdot_std];                  % Measurement Uncertainty
            
            obj.P = [1,0,0,0;                              % Initial Error Covariance Matrix
                     0, 0.02, 0, 0;
                     0, 0, 0.4,0;
                     0, 0, 0, 0.1];  
            
            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
            obj.dk_history = zeros(1,1);
            obj.P_adapt = obj.P;
            obj.wk = 0.09*[dt^2;dt;dt^2;dt];
            obj.count =0;
            obj.counter =0;
            obj.updater =0;
            obj.update1 =0;

        end
        
        function [X_pred,KF_obj1] = predict(obj)

            %PREDICT NEXT STATE (prior)
            % x = Fx
            obj.X = obj.F*obj.X +obj.wk ; 
            
            % P = FPF' + Q
            obj.P = obj.A * obj.P * obj.A.' + obj.Q;

            X_pred = obj.X;
            KF_obj1  = obj;

        end
        
        function [X_est,KF_obj2] = update(obj,z)
            %UPDATE STAGE
            obj.count = obj.count+1;
            %Adaptive Filtering
            threshold=3;
            dk = obj.X(3,1)-z(2);
            y = dk*dk;
            %S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            eps_ =log(normpdf(z(2),obj.X(3,1),obj.S(2,2)));
            disp("eps_");
            disp(eps_);

            if(abs(eps_)>threshold && obj.count>10 && obj.updater<4 && obj.update1>4)
                disp("Threshold");
                obj.updater = obj.updater+1;
                if(obj.updater==4)
                    obj.update1=0;
                end

                obj.P_adapt =obj.P;
                k_p = (obj.P_adapt * obj.H.') * obj.S^(-1);

                r_adapt = obj.R;
                r_adapt(2,2) =r_adapt(2,2)*1000;

                obj.S = obj.H * obj.P * obj.H.' + r_adapt;

                K = (obj.P * obj.H.') * obj.S^(-1);
                K(3,2) = 0;
                K(4,2) = 0;
                        
                obj.X = obj.X + K * (z-obj.H * obj.X);
                I = eye(size(obj.H,2));
        
                %UPDATE ERROR COVARIANCE MATRIX
                obj.P = (I - (k_p * obj.H)) * obj.P ;
            else
                obj.updater =0;
                obj.update1=obj.update1+1;
                obj.S = obj.H * obj.P * obj.H.' + obj.R;
                            %K = PH'inv(S)
                K = (obj.P * obj.H.') * obj.S^(-1);
                
                obj.X = obj.X + K * (z-obj.H * obj.X);
                I = eye(size(obj.H,2));
        
                %UPDATE ERROR COVARIANCE MATRIX
                obj.P = (I - (K * obj.H)) * obj.P ;
            end
 
           

          
            X_est = obj.X;
            KF_obj2 = obj;
        end

        
     end
 end
        
    


