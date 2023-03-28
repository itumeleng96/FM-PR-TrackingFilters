classdef particleFilter

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,particles,weights;
    end
    
    methods
        function obj = particleFilter(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,initialCentroid,N)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % u_x : acceleration in the x-direction
            % u_y : acceleration in the y-direction
            % x_std_meas : standard deviation of the measurement in the x-direction
            % y_std_meas : standard deviation of the measurement in the y-direction
            
            
            %Control Input Variables
            obj.U = [u_x;
                     u_y];
    
           %Initial State of the particles
            %Create Gaussian or uniformly distributed particles on Initialization
            obj.X= [0;0;0;0;0;0];
            obj.particles = createGaussianParticles(initialCentroid,[1,10],N);

            obj.weights = ones(N,1)/N;

            %State transition matrix
            obj.dt = dt;


            obj.A = [1,dt,(1/2)*dt^2, 0 ,0 ,0;
                     0, 1, dt,0,0,0;
                     0, 0, 1, 0,0,0;
                     0, 0  0, 1,dt,(1/2)*dt^2;
                     0, 0, 0, 0,1,dt ;
                     0, 0, 0, 0,0,1;];


            %The control input matrix B
            obj.B = [0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;];

            %Measure ment Mapping Matrix 
            obj.H = [1,0,0,0,0,0;
                     0,0,0,1,0,0];

            %Process Noise Covariance
            obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2, 0,0,0;
                     (dt^3)/2, dt^2, dt, 0,0,0;
                     (dt^2)/2, dt, 1, 0,0,0;
                     0, 0, 0, (dt^4)/4, (dt^3)/2, (dt^2)/2;
                     0, 0, 0, (dt^3)/2, dt^2, dt;
                     0, 0, 0, (dt^2)/2, dt,1].*std_acc^2;

            
            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.A,2));

        end
        
        function [X_pred,PF_obj] = predict(obj)
            %Calculate the predicted time state

            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            %obj.X= obj.A * obj.X + obj.B * obj.U;
            
            %particles : 6XN matrix where N is number of particles
            % Calculate the predicted time state for all particles
            X_pred = obj.A * obj.particles';
            obj.particles = X_pred';

            X_pred = mean(obj.particles,1)';

            %Get mean of particles (Prediction centroid)

            %calculate error covariance
            %P= A*P*A' + Q 
            obj.P = eye(size(obj.A,2));
            
            obj.P = (obj.A * obj.P) * obj.A.' + obj.Q;
            
            PF_obj  = obj;
        end  
        
        function [X_est,PF_obj] = update(obj,z)
            %Update stage
            obj.weights(:)= 1;
            meas_err = 1;

            %Get distance between particles and the measured values
            distance =  sqrt((obj.particles(:,1)- z(1)).^2 + (obj.particles(:,4)- z(2)).^2);
            
            %Get the shortest distance and weigh
            mean = min(distance);
            obj.weights =  obj.weights .* normpdf(distance,mean,meas_err);
            
            obj.weights = obj.weights + 1.e-300;
            obj.weights = obj.weights/(sum(obj.weights));
            
            %resample if too few effective particles,duplicate useful particles
            neff = NEFF(obj.weights);
            if neff< N/2 
                indexes = resampleSystematic(obj.weights);
                [obj.particles,obj.weights]= resampleFromIndex(obj.particles,indexes);
            end        
            
            [mean,var] = estimate(obj.particles,obj.weights);


            I = eye(size(obj.H,2));

            %Update Error Covariance matrix
            obj.P = (I - (K * obj.H)) * obj.P;
            
            X_est = obj.X;
            PF_obj = obj;
        end
     end


 end
        
    


