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
            obj.F = [1, obj.k_d*dt;
                     0, 1];



            %Process Noise Covariance Matrix For Random Acceleration
            
            obj.Q = [obj.k_d^4*(dt^4)/4,obj.k_d^2*(dt^3)/2;
                    obj.k_d^2*(dt^3)/2,dt^2]*std_acc;
            
            %Measurement Error covariance matrix
            obj.R = [r_std^2,0;
                     0,rdot_std^2];


            obj.P = eye(size(obj.F,2));

            %Constants for sigma points
            obj.alpha =0.1;
            obj.kappa =0;
            obj.beta =2;
            obj.n = 2;
            obj.lambda = obj.alpha^2*(obj.n+obj.kappa) -obj.n;
            
            [obj.Wc,obj.Wm] =obj.createWeights();
            obj.sigmaPoints = obj.createSigmaPoints(obj.X);


        end
        
        function [X_pred,UKF_obj1] = predict(obj)
            %PREDICTION STAGE
            
            %x_k = F*x_(k-1)
            %calculate sigma points for given mean and covariance
            obj.sigmaPoints = obj.createSigmaPoints(obj.X);
                      
            obj.sigmaPoints(:, 1:2) = obj.F(1:2, 1:2) * obj.sigmaPoints(:, 1:2)';

            [obj.X,obj.P] = obj.unscentedTransform();


            X_pred = obj.X;
            UKF_obj1  = obj;
        end
        
        function [Xest,KF_obj2] = update(obj,z)
            
            muZ = sum(obj.sigmaPoints .* obj.Wm,1);
            y = z-muZ';

            diff = obj.sigmaPoints -muZ;
            Pz = sum(obj.Wc .* (diff * diff'), 1) + obj.R;  

            %Kalman Gain 
            K = sum(obj.Wc .* ((obj.sigmaPoints-obj.X) *(obj.sigmaPoints - muZ)'), 3) + obj.R;  
            
            obj.X  = obj.X + K*y;
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
            
            sigmaPoints = zeros(2 * obj.n + 1, 2);
            U = chol((obj.n + obj.kappa) * obj.P, 'lower');
            
            sigmaPoints(1, :) = meanValue;
            for k = 1:obj.n
                sigmaPoints(k + 1, :) = meanValue + U(k, :);
                sigmaPoints(obj.n + k + 1, :) = meanValue - U(k, :);
            end
        end

        function [meanValue,P] = unscentedTransform(obj)
            meanValue = dot(obj.sigmaPoints,obj.Wm);

            %P=∑w(X −μ)(X −μ)T+Q
            
            P = zeros(obj.n, obj.n);
            
            % Calculate the weighted sum of outer products using matrix operations

            y_diff = obj.sigmaPoints - meanValue;
            P = P + sum(obj.Wc .* (y_diff' * y_diff), 3);  
            P = P + obj.Q;

        end


    end

 end
        
    


