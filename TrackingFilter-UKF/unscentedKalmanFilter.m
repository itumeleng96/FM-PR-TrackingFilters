classdef unscentedKalmanFilter

    properties
        dt,X,F,A,H,Q,R,P,S,
        coeff,measured_x,measured_y,std_acc,k_d,Wm,Wc,
        lambda,alpha,kappa,n,beta,sigmaPoints,count;
    end
    
    methods
        function obj = unscentedKalmanFilter(dt,std_acc,r_std,rdot_std,X_initial)
        
            %Init funtion
            %Inputs: 
            % dt : samling time
            
        
            %Initial State
            obj.X= X_initial;
                
            %Update Interval
            obj.dt = dt;

            %wave number k=-lambda=c/f
            c=299792458;
            k = -c/94e6;                                    
                    
            obj.F = [1, 0, k*dt,0;
                     0, 0, k, k*dt;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            %Process Noise Covariance Matrix For Random Acceleration
            

            obj.Q = [5,0,0,0;
                     0, 0.02, 0, 0;
                     0, 0, 0.0002,0;
                     0, 0, 0, 0.005];

            %Measurement Error covariance matrix
            obj.R = [r_std^2,0;
                     0,rdot_std^2;];


            obj.P = [5,0,0,0;                             
                     0, 0.02, 0, 0;
                     0, 0, 0.4,0;
                     0, 0, 0, 0.1];    

            obj.H = [1,0,0,0;
                     0,0,1,0;];                             % Measurement Function'

            %Constants for sigma points
            obj.alpha =1;
            obj.kappa =0;
            obj.beta =2;
            obj.n = 4;
            obj.lambda = obj.alpha^2*(obj.n+obj.kappa) -obj.n;
            
            [obj.Wc,obj.Wm] =obj.createWeights();

        end
        
        function [X_pred,UKF_obj1] = predict(obj)
            %PREDICTION STAGE
            
            %calculate sigma points for given mean and covariance
            obj.sigmaPoints = obj.createSigmaPoints(obj.X');
                      
            obj.sigmaPoints = obj.F(1:4, 1:4) * obj.sigmaPoints(:, 1:4)';

            [obj.X,obj.P] = obj.unscentedTransform();


            X_pred = obj.X';
            UKF_obj1  = obj;
        end
        
        function [Xest,KF_obj2] = update(obj,z)
            muZ = sum(obj.sigmaPoints' .* obj.Wm,1);
            Xu = sum(obj.sigmaPoints' .* obj.Wm,1);
            
            y = z-obj.H*muZ';


            Pz = obj.unscentedTransformZ(muZ);  
            Pxz=obj.unscentedTransformCross(Xu,muZ);
            obj.S =Pxz;
            
            %Kalman Gain 
            K =Pxz' * Pz^(-1) ;  
            
            obj.X  = obj.X' + K*y;
            obj.X=obj.X';
            obj.P = obj.P - K*Pz*K';
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
            U = chol((obj.n + obj.kappa) * obj.P + epsilon * eye(4));
            
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

            [num,kmax] = size(obj.sigmaPoints);
            P = 0;


                        
            for k = 1:kmax
                y_diff = obj.sigmaPoints(:, k) - meanValue';
                P = P + obj.Wc(k) * (y_diff * y_diff');
            end

            P = P + obj.Q;

        end

        function PZ = unscentedTransformZ(obj,muZ)
            
            [num,kmax] = size(obj.sigmaPoints);

            PZ = zeros(2,2);
            
            % Calculate PZ using a for loop
            for k = 1:kmax
                y_diff = obj.H*obj.sigmaPoints(:,k) - obj.H*muZ';
                PZ = PZ + obj.Wc(k) * (y_diff * y_diff');
            end
            
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
                x_diff = obj.sigmaPoints(:,k)' - stateMean;
                
                Pxz = Pxz + obj.Wc(k) * (x_diff.* y_diff);
            end
            
        end
        

    end

 end
        
    