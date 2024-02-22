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
            obj.S = obj.R;

            obj.epsDoppler =[];


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

            %Adaptive estimation of R matrix
           

            %S = H*P*H'+ R
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Covariance Scaling Kalman Step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mk = e'.[HPH' +R]^−1.ek .
            ek = z-obj.H*obj.X;
            
            ek2_doppler = ek(2)^2;
            disp("ek^2");
            disp(ek2_doppler);
            disp(obj.S(2,2));
            eps_doppler = ek2_doppler/obj.S(2,2);  

            obj.epsDoppler = [obj.epsDoppler,eps_doppler];
            
            M=5;
            if(size(obj.epsDoppler,2)>M)
                average_eps = mean(obj.epsDoppler(end-M:end-1));
                disp(average_eps);
                disp(eps_doppler);
                if(eps_doppler>2*average_eps && eps_doppler>1)
                    disp("outlier");
                    obj.R = abs([0;eps_doppler-1].*(obj.H*obj.P*obj.H') + [1;eps_doppler*500].*obj.R);
                    obj.S = obj.H * obj.P * obj.H.' + obj.R;
                    %remove outlier from average
                    obj.epsDoppler =obj.epsDoppler(end-M:end-1);
                else
                    obj.R = [500,0;0,0.01];
                    obj.S = obj.H * obj.P * obj.H.' + obj.R;
                end
            end


            %K = PH'inv(S)
            K = (obj.P * obj.H.') * obj.S^(-1);
            
            %Robust estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Represent Kalman filter as linear regression problem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate dimensions of P and R matrices
            [num_rows_P, num_cols_P] = size(obj.P);
            [num_rows_R, num_cols_R] = size(obj.R);
            

            % Create block matrix representing covariance of error'

            epsilon_covariance_range =  [obj.P(1:2,1:2) ,zeros(num_rows_P/2, num_cols_R/2); 
                                   zeros(num_rows_R/2, num_cols_P/2), obj.R(1,1)];

            epsilon_covariance_doppler = [obj.P(3:4,3:4), zeros(num_rows_P/2, num_cols_R/2); 
                                   zeros(num_rows_R/2, num_cols_P/2), obj.R(2,2)]; 
            

            Sadapt_range = chol(epsilon_covariance_range,'lower');
            Sadapt_doppler = chol(epsilon_covariance_doppler,'lower');

            Ydoppler = Sadapt_doppler\ [obj.X(3:4);z(2)];
            Yrange = Sadapt_range \ [obj.X(1:2);z(1)];
                        
            Xadapt_range = Sadapt_range\ [eye(2);[1,0;]];
            Xadapt_doppler = Sadapt_doppler\ [eye(2); [1,0;]];
        
            % Define the anonymous function for the criterion
            robustScore = HuberScore(1.5);
            robustScoreH = HuberScore(1.5);
            
            % Initial guess for the state vector
            x0_range = obj.X(1:2);
            x0_doppler = obj.X(3:4);
            
            [res_doppler, ~] = fminsearch(@(xx) m_estimate_criterion(xx, Ydoppler, Xadapt_doppler,robustScoreH), x0_doppler);

            [res_range, ~] = fminsearch(@(xx) m_estimate_criterion(xx, Yrange, Xadapt_range,robustScore), x0_range);

            
            X_res = [res_range;res_doppler];
            obj.X = X_res;

            %obj.X = obj.X + K * (z-obj.H * obj.X);
            
            I = eye(size(obj.H,2));
    
            %UPDATE ERROR COVARIANCE MATRIX
            obj.P = (I - (K * obj.H)) * obj.P ; 
            obj.R = [500,0;0,0.1];
            
            X_est = obj.X;
            KF_obj2 = obj;
        end
        
     end
 end
        
    


