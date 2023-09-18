classdef unscentedKalmanFilter

    properties
        dt,X,F,A,H,Q,R,P,S,
        coeff,measured_x,measured_y,std_acc,k_d,Wm,Wc,
        lambda,alpha,kappa,n,beta,sigmaPoints;
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
            obj.k_d = -c/94e6; 

            %State transition matrix
            obj.F = [1,obj.k_d*dt,obj.k_d*(1/2)*dt^2;
                     0, 1, dt;
                     0, 0, 1;];

            %Process Noise Covariance Matrix For Random Acceleration
            
            obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2;
                     (dt^3)/2, dt^2, dt;
                     (dt^2)/2, dt , 1]*std_acc;
            
            %Measurement Error covariance matrix
            obj.R = [r_std^2,0,0;
                     0,rdot_std^2,0;
                     0,0,0.1];


            obj.P = eye(size(obj.F,2));

            obj.S = [0,0;0,0.0];                  % System Uncertainty  

            %Constants for sigma points
            obj.alpha =1;
            obj.kappa =0;
            obj.beta =2;
            obj.n = 3;
            obj.lambda = obj.alpha^2*(obj.n+obj.kappa) -obj.n;
            
            [obj.Wc,obj.Wm] =obj.createWeights();
            obj.H = [1,0,0;0,1,0;0,0,0;];

        end
        
        function [X_pred,UKF_obj1] = predict(obj)
            %PREDICTION STAGE
            
            %calculate sigma points for given mean and covariance
            obj.sigmaPoints = obj.createSigmaPoints(obj.X');
                      
            obj.sigmaPoints = obj.F(1:3, 1:3) * obj.sigmaPoints(:, 1:3)';

            [obj.X,obj.P] = obj.unscentedTransform();


            X_pred = obj.X';
            UKF_obj1  = obj;
        end
        
        function [Xest,KF_obj2] = update(obj,z)
            z= [z;0];
            muZ = sum(obj.sigmaPoints' .* obj.Wm,1);
            Xu = sum(obj.sigmaPoints' .* obj.Wm,1);
            
            y = z-muZ';

            [~,Pz] = obj.unscentedTransformZ();  
            Pxz=obj.unscentedTransformCross(Xu,muZ);
            obj.S =Pxz;
            
            %Kalman Gain 
            K =Pxz* Pz^-1 ;  
            
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
            
            sigmaPoints = zeros(2 * obj.n + 1, 3);
            %U = chol((obj.n + obj.kappa) * obj.P);
            epsilon = 1e-6;  % Small positive perturbation
            U = chol((obj.n + obj.kappa) * obj.P + epsilon * eye(size(obj.P)));
            
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
            P = zeros(num, num);


                        
            for k = 1:kmax
                y_diff = obj.sigmaPoints(:, k)' - meanValue;
                P = P + obj.Wc(k) * (y_diff * y_diff');
            end

            P = P + obj.Q;

        end

        function [meanValueZ,PZ] = unscentedTransformZ(obj)
            meanValueZ = sum(obj.sigmaPoints' .* obj.Wm,1);
            [num,kmax] = size(obj.sigmaPoints);

            PZ = zeros(num,num);
            
            % Calculate PZ using a for loop
            for k = 1:kmax
                y_diff = obj.sigmaPoints(:,k)' - meanValueZ;
                PZ = PZ + obj.Wc(k) * (y_diff * y_diff');
            end
            
            % Add the measurement noise covariance matrix R
            PZ = PZ + obj.R;

        end

        function [Pxz] = unscentedTransformCross(obj,stateMean,measMean)

            %P=∑w(X −μ)(X −μ)T+Q
            
            [num,kmax] = size(obj.sigmaPoints);

            Pxz = zeros(num,num);
            
            % Calculate PZ using a for loop
            for k = 1:kmax
                y_diff = obj.sigmaPoints(:,k)' - measMean;
                x_diff = obj.sigmaPoints(:,k)' - stateMean;
                
                Pxz = Pxz + obj.Wc(k) * (y_diff * x_diff');
            end
            
        end
        

    end

 end
        
    


