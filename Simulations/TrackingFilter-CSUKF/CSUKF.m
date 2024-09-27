classdef CSUKF

    properties
        dt,X,F,A,H,Q,R,P,Pk,S,
        coeff,measured_x,measured_y,std_acc,k_d,Wm,Wc,
        lambda,alpha,kappa,n,beta,sigmaPoints,wk,count,updater,update1,epsDoppler,epsRange;
    end
    
    methods
        function obj = CSUKF(dt,std_acc,r_std,rdot_std,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : samling time
            
        
            %Initial State
            obj.X= X_initial;
                
            %Update Interval
            obj.dt = dt;

           obj.F = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];                            % Measurement Function

            

            obj.Q = std_acc*[(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            %Measurement Error covariance matrix
            obj.R = [r_std,0;
                     0,rdot_std;];

            obj.P = [5,0,0,0;                                      % Initial Error Covariance Matrix
                     0, 2.5, 0, 0;
                     0, 0, 1,0;
                     0, 0, 0, 0.5];    

            obj.H = [1,0,0,0;
                     0,0,1,0;];                                         % Measurement Function'

            %Tuning parameters For UKF
            obj.alpha =0.01;                 %Determines the spread of the sigma points around the mean : usually + value : a=0.0001
            obj.kappa =0;                   %Secondary scaling parameter usually  set to zero 
            obj.beta =2;                    %Is used to incorporate prior knowledge of the distribution of the input random variable (Gaussian : beta=2)
            obj.n = 4;                      %Is the number of dimensions 

            obj.lambda = obj.alpha^2*(obj.n+obj.kappa) -obj.n;
            
            [obj.Wc,obj.Wm] =obj.createWeights();
            obj.wk = std_acc*[dt^2;dt;dt^2;dt];

            obj.epsDoppler =[];
            obj.epsRange =[];


        end
        
        function [X_pred,UKF_obj1] = predict(obj)
            %PREDICTION STAGE
            %noise = mvnrnd([0, 0, 0, 0], obj.Q,9);

            %calculate sigma points for given mean and covariance
            obj.sigmaPoints = obj.createSigmaPoints(obj.X');
                      
            %obj.sigmaPoints = obj.F(1:4, 1:4) * obj.sigmaPoints(:, 1:4)'+ noise(:, 1:4)';
            obj.sigmaPoints = obj.F(1:4, 1:4) * obj.sigmaPoints(:, 1:4)';

            [obj.X,obj.Pk] = obj.unscentedTransform();

            %to get S on predict  
            muZ = sum(obj.sigmaPoints' .* obj.Wm,1);
            Pz = obj.unscentedTransformZ(muZ);  
            obj.S = Pz;

            X_pred = obj.X' +obj.wk;
            UKF_obj1  = obj;
        end
        
        function [Xest,KF_obj2] = update(obj,z)

            muZ = sum(obj.sigmaPoints' .* obj.Wm,1);
            Xu = sum(obj.sigmaPoints' .* obj.Wm,1);
            
            y = z-obj.H*muZ';

            Pz = obj.unscentedTransformZ(muZ);  
            Pxz=obj.unscentedTransformCross(Xu,muZ);
            
            obj.S = Pz;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Covariance Scaling Unscented Kalman Filter Step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            eps_doppler = (y(2)^2)/obj.S(2,2);  
            eps_range = (y(1)^2)/obj.S(1,1);   

            obj.epsDoppler = [obj.epsDoppler,eps_doppler];
            obj.epsRange = [obj.epsRange,eps_range];
                       
            M=6; %Number of samples to average
            r_adapt= obj.R;

            if(size(obj.epsDoppler,2)>M && eps_doppler>1 && eps_doppler>mean(obj.epsDoppler(end-M:end-1)))                
                r_adapt(2,2) =(eps_doppler-1)*obj.P(3,3) + eps_doppler*r_adapt(2,2)*1000;
                %remove outlier from average
                obj.epsDoppler =obj.epsDoppler(end-M:end-1);
            end
            if(size(obj.epsRange,2)>M && eps_range>2 && eps_range>mean(obj.epsRange(end-M:end-1)))                
                r_adapt(1,1) =(eps_range-1)*obj.P(1,1) + eps_range*r_adapt(1,1);
                %remove outlier from average
                obj.epsRange =obj.epsRange(end-M:end-1);
            end

            obj.R = r_adapt;
            Pz = obj.unscentedTransformZ(muZ);  

            %{
            [num_rows_P, num_cols_P] = size(obj.P);
            [num_rows_R, num_cols_R] = size(obj.R);
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Huber's M estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Create block matrix representing covariance of error'

            epsilon_covariance_range =  [obj.P(1:2,1:2) ,zeros(num_rows_P/2, num_cols_R/2); 
                                   zeros(num_rows_R/2, num_cols_P/2), obj.R(1,1)];

            epsilon_covariance_doppler = [obj.P(3:4,3:4), zeros(num_rows_P/2, num_cols_R/2); 
                                   zeros(num_rows_R/2, num_cols_P/2), obj.R(2,2)]; 
            

            Sadapt_range = chol(epsilon_covariance_range,'lower');
            Sadapt_doppler = chol(epsilon_covariance_doppler,'lower');

            Ydoppler = Sadapt_doppler\ [obj.X(3:4)';z(2)];
            Yrange = Sadapt_range \ [obj.X(1:2)';z(1)];
                        
            Xadapt_range = Sadapt_range\ [eye(2);[1,0;]];
            Xadapt_doppler = Sadapt_doppler\ [eye(2); [1,0;]];
        
            % Define the anonymous function for the criterion
            robustScore = HuberScore(1.5);
            robustScoreH = HuberScore(1.5);
            
            % Initial guess for the state vector-
            x0_range = obj.X(1:2)';
            x0_doppler = obj.X(3:4)';
            
            [res_doppler, ~] = fminsearch(@(xx) m_estimate_criterion(xx, Ydoppler, Xadapt_doppler,robustScoreH), x0_doppler);

            [res_range, ~] = fminsearch(@(xx) m_estimate_criterion(xx, Yrange, Xadapt_range,robustScore), x0_range);

            
            X_res = [res_range;res_doppler];
            obj.X = X_res;
            
            %}
            K =Pxz * Pz^(-1) ;  
            obj.X  = obj.X' + K*y;

            obj.P = obj.Pk - K*Pz*K';
            obj.R = [500,0;0,0.1];


            obj.X=obj.X';

            %Investigate why does Doppler estimation jump around when
            %measurement unlikely
            Xest=0;
            KF_obj2 = obj;
        end
        
        function [Wc,Wm] = createWeights(obj)
            %Compute Weights according to Van Der Merwe Implementation
            
            Wc = ones(1, 2 * obj.n + 1) * (1 / (2 * (obj.n + obj.lambda)));
            Wm = ones(1, 2 * obj.n + 1) * (1 / (2 * (obj.n + obj.lambda)));

            Wm(1) = obj.lambda / (obj.n + obj.lambda);
            Wc(1) = obj.lambda / (obj.n + obj.lambda) + (1 - obj.alpha^2 + obj.beta);

            Wm=Wm';
            Wc=Wc';
            

        end

        function [sigmaPoints]  = createSigmaPoints(obj,meanValue)
            %Compute the Sigma Points
            obj.n = obj.n;
            obj.P = obj.P;
            
            sigmaPoints = zeros(2 * obj.n + 1, 4);
            %U = chol((obj.n + obj.kappa) * obj.P);
            epsilon = 1e-6;  % Small positive perturbation
            U = chol((obj.n + obj.lambda) * obj.P + epsilon * eye(4));
            
            sigmaPoints(1, :) = meanValue;

            for k = 1:obj.n
                sigmaPoints(k + 1, :) = meanValue' + U(k, :);
                sigmaPoints(obj.n + k + 1, :) = meanValue' - U(k, :);
            end
        end

        function [meanValue,P] = unscentedTransform(obj)
            meanValue = sum(obj.sigmaPoints' .* obj.Wm, 1);

            %P=∑w(X −μ)(X −μ)T+Q
                        
            % Calculate the weighted sum of outer products using matrix operations

            y_diff = obj.sigmaPoints - meanValue';
            P = sum(obj.Wc' .* y_diff * y_diff',3);
            P = P + obj.Q;

        end

        function PZ = unscentedTransformZ(obj,muZ)
                       
            y_diff = obj.H*obj.sigmaPoints - obj.H*muZ';
            PZ = sum(obj.Wc' .* y_diff * y_diff',3);

            % Add the measurement noise covariance matrix R
            PZ = PZ + obj.R;

        end

        function [Pxz] = unscentedTransformCross(obj,stateMean,measMean)

            %P=∑w(X −μ)(X −μ)T+Q
            
            [num,kmax] = size(obj.sigmaPoints);

            Pxz = 0;
            
            % Calculate PZ using a for loop
            for k = 1:kmax
                y_diff = obj.H*obj.sigmaPoints(:,k) - obj.H*measMean';
                x_diff = obj.sigmaPoints(:,k) - stateMean';
                
                Pxz = Pxz + obj.Wc(k) * (x_diff.* y_diff');
            end
            
        end
        

    end

 end
        
    