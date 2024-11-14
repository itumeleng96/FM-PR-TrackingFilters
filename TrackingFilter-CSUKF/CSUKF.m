classdef CSUKF

    properties
        dt,X,F,A,H,Q,R,P,Pk,S,
        coeff,measured_x,measured_y,std_acc,k_d,Wm,Wc,
        lambda,alpha,kappa,n,beta,sigmaPoints,wk,count,updater,update1,epsDoppler,epsRange;
    end
    
    methods
        function obj = CSUKF(dt,std_acc,r_std,rdot_std,X_initial,alpha,kappa,beta)
        
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

            

            obj.Q = [std_acc(1)*(dt^4)/4,std_acc(1)*(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     std_acc(1)*(dt^3)/2, std_acc(1)*dt^2, 0, 0;
                     0, 0, std_acc(2)*(dt^4)/4,std_acc(2)*(dt^3)/2;
                     0, 0, std_acc(2)*(dt^3)/2, std_acc(2)*dt^2];


            %Measurement Error covariance matrix
            obj.R = [r_std,0;
                     0,rdot_std;];

            obj.P = [5,0,0,0;                              % Initial Error Covariance Matrix
                     0, 1, 0, 0;
                     0, 0, 2,0;
                     0, 0, 0, 1];  
            

            obj.H = [1,0,0,0;
                     0,0,1,0;];                                         % Measurement Function'

            %Tuning parameters For UKF
              %Tuning parameters For UKF
            obj.alpha =alpha;                %Determines the spread of the sigma points around the mean : usually + value : a=0.0001
            obj.kappa =kappa;                %Secondary scaling parameter usually  set to zero 
            obj.beta =beta;                  %Is used to incorporate prior knowledge of the distribution of the input random variable (Gaussian : beta=2)
            obj.n = 4;                      %Is the number of dimensions 

            obj.lambda = obj.alpha^2*(obj.n+obj.kappa) -obj.n;
            
            [obj.Wc,obj.Wm] =obj.createWeights();
            obj.wk = [std_acc(1)*dt^2;std_acc(1)*dt;std_acc(2)*dt^2;std_acc(2)*dt];

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
        
        function [Xest, KF_obj2] = update(obj, z)
            %UPDATE STAGE
        
            % Store the original R matrix before any modifications
            original_R = obj.R;
        
            % Calculate the mean of sigma points
            muZ = sum(obj.sigmaPoints' .* obj.Wm, 1);
            Xu = sum(obj.sigmaPoints' .* obj.Wm, 1);
            
            % Measurement residual
            y = z - obj.H * muZ';
        
            % Unscented transform for measurement covariance and cross-covariance
            Pz = obj.unscentedTransformZ(muZ);
            Pxz = obj.unscentedTransformCross(Xu, muZ);
        
            obj.S = Pz;
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Covariance Scaling Unscented Kalman Filter Step 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Separate range and Doppler residuals
            eps_range = (y(1)^2) / obj.S(1,1);  % Range residual (x component)
            eps_doppler = (y(2)^2) / obj.S(2,2);  % Doppler residual (y component)
        
            % Store residuals in moving window arrays
            obj.epsRange = [obj.epsRange, eps_range];
            obj.epsDoppler = [obj.epsDoppler, eps_doppler];
        
            % Number of samples to average for outlier detection
            M = 6;
            r_adapt = obj.R;
            alphaFactor = 0.2;
        
            % Outlier rejection for Doppler (y component)
            if (size(obj.epsDoppler, 2) > M)
                residual_mean_doppler = mean(obj.epsDoppler(end-M:end-1));
                residual_std_doppler = std(obj.epsDoppler(end-M:end-1));
        
                if eps_doppler > 1 && eps_doppler > (residual_mean_doppler + residual_std_doppler)
                    r_adapt(2,2) = alphaFactor * r_adapt(2,2) + (1 - alphaFactor) * (eps_doppler + obj.P(3,3));
                    % Remove old values from moving window
                    obj.epsDoppler = obj.epsDoppler(end-M:end-1);
                end
            end
        
            % Outlier rejection for Range (x component)
            if (size(obj.epsRange, 2) > M)
                residual_mean_range = mean(obj.epsRange(end-M:end-1));
                residual_std_range = std(obj.epsRange(end-M:end-1));
        
                if eps_range > 2 && eps_range > (residual_mean_range + residual_std_range)
                    r_adapt(1,1) = alphaFactor * r_adapt(1,1) + (1 - alphaFactor) * (eps_range + obj.P(1,1));
                    % Remove old values from moving window
                    obj.epsRange = obj.epsRange(end-M:end-1);
                end
            end
        
            % Update the R matrix with the adapted values
            obj.R = r_adapt;
        
            % Recompute Pz after updating R
            Pz = obj.unscentedTransformZ(muZ);
        
            % Calculate Kalman gain
            K =Pxz * Pz^(-1) ;  
        
            % Update state estimate
            obj.X = obj.X' + K*y;
        
            % Update error covariance
            obj.P = obj.Pk - K*Pz*K';
        
            % Restore the original R matrix
            obj.R = original_R;
        
            % Return updated state estimate and object
            obj.X=obj.X';  
            Xest = obj.X;
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
        
    