classdef particleFilter

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,particles,weights,N;
    end
    
    methods
        function obj = particleFilter(dt,u_x,u_y,std_acc,initialCentroid,N)
        
            %Init funtion
            %Inputs: 
            % dt : sampling time
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
            obj.particles = obj.createUniformParticles([0,1e-4],[0,150],N);

            obj.weights = ones(N,1)/N;

            %State transition matrix
            obj.dt = dt;


            obj.A = [1,dt,(1/2)*dt^2, 0 ,0 ,0;
                     0, 1, dt,0,0,0;
                     0, 0, 1, 0,0,0;
                     0, 0  0, 1,dt,(1/2)*dt^2;
                     0, 0, 0, 0,1,dt ;
                     0, 0, 0, 0,0,1;];


            %Process Noise Covariance
            obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2, 0,0,0;
                     (dt^3)/2, dt^2, dt, 0,0,0;
                     (dt^2)/2, dt, 1, 0,0,0;
                     0, 0, 0, (dt^4)/4, (dt^3)/2, (dt^2)/2;
                     0, 0, 0, (dt^3)/2, dt^2, dt;
                     0, 0, 0, (dt^2)/2, dt,1];

        
            obj.N = N;

        end
        
        function [X_pred,PF_obj] = predict(obj)
            %Calculate the predicted time state

            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            %obj.X= obj.A * obj.X + obj.B * obj.U;
            
            %particles : 6XN matrix where N is number of particles
            % Predict new particle states by adding Gaussian noise to each particle
            noise = mvnrnd(zeros(size(obj.Q,1),1), obj.Q, obj.N);
            obj.particles = obj.A * obj.particles' + noise';
            obj.particles = obj.particles';
            
            X_pred = mean(double(obj.particles),1)';
            disp("Prediction");
            disp(X_pred);
            %Get mean of particles (Prediction centroid)

            PF_obj  = obj;
        end  
        
        function [X_est,PF_obj] = update(obj,z)
            %Update stage
            obj.weights(:)= 1;
            meas_err = 0.01;

            %Get distance between particles and the measured values
            distance =  sqrt((obj.particles(:,1)- z(1)).^2 + (obj.particles(:,4)- z(2)).^2);
            
            %Get the shortest distance and weigh
            mean = min(distance);
            obj.weights =  obj.weights .* normpdf(distance,mean,meas_err);
            
            obj.weights = obj.weights + 1.e-300;
            obj.weights = obj.weights/(sum(obj.weights));
            
            %resample if too few effective particles,duplicate useful particles
            neff = obj.NEFF(obj.weights);
            if neff< obj.N/2 
                indexes = obj.resampleSystematic(obj.weights);
                [obj.particles,obj.weights]= obj.resampleFromIndex(obj.particles,indexes);
            end        
            
            [mean,var] = obj.estimate(obj.particles,obj.weights);
            disp(mean);
            %Update Error Covariance matrix            
            X_est = mean;
            PF_obj = obj;
        end
    end
    methods(Static)
        function [ indx ] = resampleSystematic( w )
            N = length(w);
            Q = cumsum(w);
            T = linspace(0,1-1/N,N) + rand(1)/N;
            T(N+1) = 1;
            i=1;
            j=1;
            while (i<=N)
                if (T(i)<Q(j))
                    indx(i)=j;
                    i=i+1;
                else
                    j=j+1;        
                end
            end
        end
        function [particles] = createUniformParticles(x_range,y_range,N)
            %Create a uniform Distribution of particles over a region
            % N : number of particles
                particles = zeros(N,6);
                particles(:,1) = unifrnd(x_range(1),x_range(2),[N 1]);
                particles(:,4) = unifrnd(y_range(1),y_range(2),[N 1]);
        end
            
            
        function [particles] = createGaussianParticles(mean,std,N)
            %Create a Gaussian Distribution of particles over a region
            % N : number of particles
                particles = zeros(N,6);
                particles(:,1) = mean(1) + (randn(N,1))*std(1) ; 
                particles(:,4) = mean(2) + (randn(N,1))*std(2) ;
                    
        end
        
        function [particles,weights] = resampleFromIndex(particles,indexes)
            %RESAMPLEFROMINDEX Summary of this function goes here
            %   Detailed explanation goes here
                particles(:,1) = particles(indexes,1);
                particles(:,2) = particles(indexes,2);
                particles(:,3) = particles(indexes,3);
                particles(:,4) = particles(indexes,4);
                particles(:,5) = particles(indexes,5);
                particles(:,6) = particles(indexes,6);
                
                N = size(particles,1);
                weights = zeros(N,1);
                weights(:) = 1.0/size(weights,1);
        end
            
        function [mean, var] = estimate(particles, weights)
            %ESTIMATE Summary of this function goes here
            %   Detailed explanation goes here
            
                % Initialize mean and var to zero vectors of size 1x2
                mean = [0,0];
                var = [0,0];
                
                % Calculate the mean of particles using only columns 1 and 4
                mean(1) = sum(particles(:, 1).*weights)/sum(weights);
                mean(2) = sum(particles(:, 4).*weights)/sum(weights);
                
                % Calculate the variance of particles using only columns 1 and 4
                var_particles = zeros(size(particles, 1),2);
                var_particles(:, 1) = (particles(:, 1) - mean(1)).^2;
                var_particles(:, 2) = (particles(:, 4) - mean(2)).^2;
                
                var(1) = sum(var_particles(:, 1).*weights)/sum(weights);
                var(2) = sum(var_particles(:, 2).*weights)/sum(weights);
        end
            

        function [neff] = NEFF(weights)
            %NEFF Summary of this function goes here
            %   Detailed explanation goes here
                neff = 1./ sum(weights.^2) ;
        end
            
     end


 end
        
    


