%CSUKF  Covariance-Scaling Unscented Kalman Filter for bistatic tracking.
%
%   State vector:   x = [r, rdot, f, fdot]'
%   Measurement:    z = [r, f]'
%   Motion model:   Nearly Constant Velocity (NCV) with sampling interval dt
%
%   Propagates the state through 2n+1 = 9 sigma points using the unscented
%   transform. Adaptive covariance scaling (Akhlaghi 2018) is applied to R
%   (always) and to Q (when obj.cs_adapt_q = true) on every accepted
%   measurement, with forgetting factor obj.cs_alpha (default 0.98).
%
%   Construction:
%     f = CSUKF(dt, std_acc, r_std, rdot_std, X_initial, alpha, kappa, beta)
%       dt         - scan interval (s)
%       std_acc    - [ax, ay] process-noise std for range and Doppler axes
%       r_std      - measurement std, range channel
%       rdot_std   - measurement std, Doppler channel
%       X_initial  - initial state vector [r0; rdot0; f0; fdot0]
%       alpha, kappa, beta  - UKF sigma-point spread parameters
%                              (typical: alpha=0.5, kappa=0, beta=2)
%
%   Public methods:
%     [xhat, obj] = predict(obj)      sigma-point prediction
%     [xhat, obj] = update(obj, z)    measurement update on z = [r; f]
%
%   Adaptive tuning knobs (set after construction):
%     obj.cs_alpha    - forgetting factor in [0,1]           (default 0.98)
%     obj.cs_adapt_q  - true = also adapt Q, false = R only  (default false)
%
%   Reference:
%     Julier & Uhlmann (1997) for UKF; Akhlaghi et al. (2018) for covariance
%     adaptation. See CSKF for the R/Q update formulae.
classdef CSUKF

    properties
        dt, X, F, H, Q, Q_nominal, R, R_nominal, P, Pk, S, std_acc, ...
        Wm, Wc, lambda, alpha, kappa, n, beta, sigmaPoints, wk, ...
        cs_alpha, cs_adapt_q;
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
            obj.R_nominal = obj.R;                 % Nominal R for fading-memory decay

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

            % Akhlaghi (2018) covariance-scaling defaults (settable post-construction).
            obj.cs_alpha   = 0.95;   % forgetting factor (both R and Q)
            obj.cs_adapt_q = true;   % master switch for Q adaptation
            obj.Q_nominal  = obj.Q;  % nominal Q retained as PSD floor


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
            % Akhlaghi (2018) verbatim continuous R + Q. Outlier rejection is
            % handled upstream at the association step (GNN gate in runner).
            %   R <- alpha R + (1-alpha)( eps eps' + H P^- H' )       Eq. 11
            %   Q <- alpha Q + (1-alpha)( K d d' K' )                 Eq. 15

            alphaFactor = obj.cs_alpha;

            % Weighted mean of sigma points (Van der Merwe UKF).
            muZ = sum(obj.sigmaPoints' .* obj.Wm, 1);
            Xu  = sum(obj.sigmaPoints' .* obj.Wm, 1);

            % Innovation and unscented-transform covariances (pre-update).
            y     = z - obj.H * muZ';                   % d_k
            Pz    = obj.unscentedTransformZ(muZ);       % Pz = H Sigma H' + R
            Pxz   = obj.unscentedTransformCross(Xu, muZ);
            obj.S = Pz;

            P_pred = obj.Pk;                            % predicted P for R update

            % Standard UKF update.
            K     = Pxz * Pz^(-1);
            obj.X = obj.X' + K * y;
            obj.P = obj.Pk - K * Pz * K';

            % Akhlaghi Eq. 11 — continuous R update.
            eps_post = z - obj.H * obj.X;
            HPH      = obj.H * P_pred * obj.H.';
            obj.R    = alphaFactor * obj.R + (1 - alphaFactor) * ...
                       (eps_post * eps_post.' + HPH);

            % Akhlaghi Eq. 15 — continuous Q update (per block).
            if obj.cs_adapt_q
                KyyK = K * (y * y.') * K.';
                obj.Q(1:2,1:2) = alphaFactor * obj.Q(1:2,1:2) + ...
                                 (1 - alphaFactor) * KyyK(1:2,1:2);
                obj.Q(3:4,3:4) = alphaFactor * obj.Q(3:4,3:4) + ...
                                 (1 - alphaFactor) * KyyK(3:4,3:4);
                obj.Q(1:2,1:2) = max(obj.Q(1:2,1:2), obj.Q_nominal(1:2,1:2));
                obj.Q(3:4,3:4) = max(obj.Q(3:4,3:4), obj.Q_nominal(3:4,3:4));
            end

            obj.X   = obj.X';
            Xest    = obj.X;
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
        
    