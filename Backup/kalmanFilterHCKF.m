classdef kalmanFilter

    properties
        dt,U,X,F,A,H,Q,R,P,S,coeff,measured_x,measured_y,std_acc,wk,epsDoppler;
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
            k = -c/94e6;                                   % Wave number k=-lambda=c/f

                    
            obj.F = [1, 0, k*dt,0;
                     0, 0, k, k*dt;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];                            % Measurement Function

            


            obj.Q = [(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            obj.R = [r_std,0;
                     0,rdot_std];                  % Measurement Uncertainty
            
            obj.P = [500,0,0,0;                              % Initial Error Covariance Matrix
                     0, 10, 0, 0;
                     0, 0, 1,0;
                     0, 0, 0, 0.1];  

            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];

            obj.wk = std_acc*[dt^2;dt;dt^2;dt];

            obj.epsDoppler =[];

            obj.S = obj.R;

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

            %S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Covariance Scaling Kalman Step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mk = e'.[HPH' +R]^−1.ek .
            ek = z-obj.H*obj.X;
            
            ek2_doppler = ek(2)^2;
            eps_doppler = ek2_doppler/obj.S(2,2);  

            obj.epsDoppler = [obj.epsDoppler,eps_doppler];
            
            M=5;

            if(size(obj.epsDoppler,2)>M)
                average_eps = mean(obj.epsDoppler(end-M:end-1));
                disp(average_eps);
                if(eps_doppler>2*average_eps)
                    disp("Outlier");
                    disp(eps_doppler);
                    obj.R = [0;eps_doppler-1].*(obj.H*obj.P*obj.H') + [1;eps_doppler].*obj.R;
                    obj.R = abs(obj.R);
                    obj.S = obj.H * obj.P * obj.H.' + obj.R;
                    %remove 
                    obj.epsDoppler =obj.epsDoppler(end-M:end-1);
                else
                    obj.R = [500,0;0,0.1];
                end
            end

            %Scaling of the Measurement Covariance
            %if(eps_doppler>1)
            %    obj.R = [0;eps_doppler-1].*(obj.H*obj.P*obj.H') + [1;eps_doppler].*obj.R;
            %    obj.S = obj.H * obj.P * obj.H.' + obj.R;

            %           else
            %              obj.R = [500,0;0,0.1];
            %          end

            disp(obj.R);
           
            %K = PH'inv(S)
            K = (obj.P * obj.H.') * obj.S^(-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Huber-Kalman step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % KF as linear regression problem

            [num_rows_P, num_cols_P] = size(obj.P);
            [num_rows_R, num_cols_R] = size(obj.R);

            epsilon_covariance = [obj.P, zeros(num_rows_P, num_cols_R); 
                                   zeros(num_rows_R, num_cols_P), obj.R];
    
            Sadapt = chol(epsilon_covariance);
            Sinv = Sadapt \ eye(size(Sadapt));
            Y = Sinv * [obj.X;z];
            Xadapt = Sinv * [eye(size(obj.X,1)); obj.H];
            robustScore = HuberScore(1.5);
            
            % Initial guess for the state vector
            x0 = obj.X;

            % Perform optimization using fminsearch
            [x_min, ~] = fminsearch(@(xx) m_estimate_criterion(xx, Y, Xadapt,robustScore), x0);
            
            obj.X = x_min;

            %obj.X = obj.X + K * (z-obj.H * obj.X);
            
            I = eye(size(obj.H,2));
    
            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P ; 
           
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


