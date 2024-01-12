classdef particleFilter

    properties
        dt,             %%Sampling Time
        F,              %The state transition matrix
        Q,              %The Process Noise Covariance Matrix
        particles,      %Matrix containing the current State of the Particles
        weights,        %A vector containig the current weights of the particles
        N,              %Number of particles to use 
        scaling_factor; % 
        std_meas;
        S;
        count;
        updater;
        update1;

    end
    
    methods
        function obj = particleFilter(dt,std_acc,std_meas,initialCentroid,N)
        
            %Init funtion
            obj.N = N;
            obj.std_meas=std_meas;
           
            %Create  Gaussian  distributed particles on Initialization
            obj.particles = obj.createGaussianParticles(initialCentroid,[1000,10],N);
            
            %Set Equal weights
            obj.weights = ones(N,1)/N;
            obj.dt = dt;

            %wave number c/f
            k = -299792458/94e6; 

            obj.F = [1, 0, k*dt,0;
                     0, 0, k, k*dt;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
           
            

            obj.Q = [(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];
            
            obj.count =0;
            obj.updater =0;
            obj.update1 =0;




        end
        
        function [X_pred, PF_obj] = predict(obj)

            % Generate random Gaussian noise with zero mean and covariance matrix Q
            noise = mvnrnd([0, 0, 0, 0], obj.Q, obj.N);


            % Add process noise to particle states
            obj.particles(:, 1:4) = (obj.F(1:4, 1:4) * obj.particles(:, 1:4)' + noise(:, 1:4)')';
                    
            X_pred = mean(obj.particles, 1)';

            PF_obj = obj;
        end
        
        function [X_est, PF_obj] = update(obj, z)
        
        
            % Calculate the particle likelihoods based on a Gaussian PDF
            % p(z t​∣x t(i))=p(zx∣x t(i),σx)⋅p(z y∣y t(i),σy)
            %Adaptive filtering

            % Calculate the cross-covariance elements (K_ij) using the weighted particles
            % S = HPH' + R;
            threshold=2;

            %max number of adaptive filtering
            max_adapt=3;

            %%To be replaced with another scheme to detect steady state
            obj.count = obj.count+1;    

            
            [meanValueS,~] = obj.estimate(obj.particles, obj.weights);
            slikelihood= sum(obj.weights .* (obj.particles(:, [1 3]) - meanValueS) .* (z - meanValueS')', 3) + obj.std_meas;
            obj.S = [mean(slikelihood(:,1)),0;0,mean(slikelihood(:,2));];
            eps_ = log(normpdf(z(2),meanValueS(2),obj.S(2,2)));

            if(abs(eps_)>threshold && obj.count>10 && obj.updater<max_adapt && obj.update1>max_adapt)
                obj.updater = obj.updater+1;
                if(obj.updater==max_adapt)
                    obj.update1=0;
                end

                diffs = (obj.particles(:, [1 3])' - z)';
                likelihood_x = exp(-0.5 * (diffs(:, 1).^2) / obj.std_meas(1)^2*50);
                likelihood_y = exp(-0.5 * (diffs(:, 2).^2) / obj.std_meas(2)^2);
                likelihood = likelihood_x .* likelihood_y;
                weights_adapt =obj.weights .*likelihood;


                [meanValue,~] = obj.estimate(obj.particles, weights_adapt);
                
            else
                obj.updater =0;
                obj.update1=obj.update1+1;
                    
                diffs = (obj.particles(:, [1 3])' - z)';
                likelihood_x = exp(-0.5 * (diffs(:, 1).^2) / obj.std_meas(1)^2);
                likelihood_y = exp(-0.5 * (diffs(:, 2).^2) / obj.std_meas(2)^2);
                likelihood = likelihood_x .* likelihood_y;
                
                % Normalize the likelihood
                likelihood = likelihood / sum(likelihood);
                
                
                % Update the particle weights
                obj.weights = obj.weights .* likelihood;
                obj.weights = obj.weights + 1.e-300;  % Add small constant to avoid zero weights
                obj.weights = obj.weights / sum(obj.weights);
            
                % Resample if too few effective particles
                neff = obj.NEFF(obj.weights);
                if neff < obj.N / 2
                    indexes = obj.resampleSystematic(obj.weights);
                    [obj.particles, obj.weights] = obj.resampleFromIndex(obj.particles, indexes);
                end

                [meanValue,~] = obj.estimate(obj.particles, obj.weights);


            end

            X_est = meanValue;            
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

                
        function [particles] = createGaussianParticles(mean,std,N)
            %Create a Gaussian Distribution of particles over a region
            % N : number of particles
                particles = zeros(N,4);
                particles(:,1) = mean(1) + (randn(N,1))*std(1) ; 
                particles(:,3) = mean(3) + (randn(N,1))*std(2) ;
                    
        end
        
        function [particles,weights] = resampleFromIndex(particles,indexes)
            %RESAMPLEFROMINDEX Summary of this function goes here
            %   Detailed explanation goes here
                particles(:,1) = particles(indexes,1);
                particles(:,2) = particles(indexes,2);
                particles(:,3) = particles(indexes,3);
                particles(:,4) = particles(indexes,4);
                
                
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
                
                % Calculate the mean of particles using only columns 1 and 2
                mean(1) = sum(particles(:, 1).*weights)/sum(weights);
                mean(2) = sum(particles(:, 3).*weights)/sum(weights);
                
                % Calculate the variance of particles using only columns 1 and 2
                var_particles = zeros(size(particles, 1),2);
                var_particles(:, 1) = (particles(:, 1) - mean(1)).^2;
                var_particles(:, 2) = (particles(:, 2) - mean(2)).^2;
                
                var(1) = sum(var_particles(:, 1).*weights)/sum(weights);
                var(2) = sum(var_particles(:, 2).*weights)/sum(weights);
        end
            
        function [neff] = NEFF(weights)
            %NEFF : Effective particle sample size
            % To measure degeneracy of the  particles
            neff = 1./ sum(weights.^2) ;
        end   

        function [likelihood] = calculateLikelihood(mean,distances,meas_err)
            %Calculate the Likelihood using a Normal Distribution
            likelihood = normpdf(distances,mean,meas_err);
        end

     end
 end