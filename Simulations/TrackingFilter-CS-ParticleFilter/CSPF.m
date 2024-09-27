classdef CSPF

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
        epsDoppler;C

    end
    
    methods
        function obj = CSPF(dt,std_acc,std_meas,initialCentroid,N)
        
            %Init funtion
            obj.N = N;
            obj.std_meas=std_meas;
           
            %Create  Gaussian  distributed particles on Initialization
            obj.particles = obj.createGaussianParticles(initialCentroid,[5,5],N);
            
            %Set Equal weights
            obj.weights = ones(N,1)/N;
            obj.dt = dt;

            %wave number c/f
            obj.F = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
           
            

            obj.Q = std_acc*[(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            
            obj.count =0;
            obj.updater =0;
            obj.update1 =0;

            obj.epsDoppler =[];



        end
        
        function [X_pred, PF_obj] = predict(obj)

            % Generate random Gaussian noise with zero mean and covariance matrix Q
            noise = randn(obj.N, 4) * chol(obj.Q);


            % Add process noise to particle states
            obj.particles(:, 1:4) = (obj.F(1:4, 1:4) * obj.particles(:, 1:4)' + noise(:, 1:4)')';
                    
            X_pred = mean(obj.particles, 1)';
            

            [meanValueS,~] = obj.estimate(obj.particles, obj.weights);
            covariance_matrix = (obj.weights .* (obj.particles(:, [1 3]) - meanValueS))' * (obj.weights.*(obj.particles(:, [1 3]) - meanValueS));
            
            slikelihood= covariance_matrix+ obj.std_meas;
            obj.S = [mean(slikelihood(:,1)),0;0,mean(slikelihood(:,2));];
            PF_obj = obj;


        end
        
        function [X_est, PF_obj] = update(obj, z)
        
        
            % Calculate the particle likelihoods based on a Gaussian PDF
            % p(z t​∣x t(i))=p(zx∣x t(i),σx)⋅p(z y∣y t(i),σy)
            %Adaptive filtering

            % Calculate the cross-covariance elements (K_ij) using the weighted particles
            % S = HPH' + R;
            
            [meanValueS,~] = obj.estimate(obj.particles, obj.weights);
            covariance_matrix = (obj.weights .* (obj.particles(:, [1 3]) - meanValueS))' * (obj.weights.*(obj.particles(:, [1 3]) - meanValueS));
            
            slikelihood= covariance_matrix+ obj.std_meas;
            obj.S = [mean(slikelihood(:,1)),0;0,mean(slikelihood(:,2));];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Detect outliers and reject them
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ek = z-meanValueS';
            ek2_doppler = ek(2)^2;

            eps_doppler = ek2_doppler/obj.S(2,2);  

            obj.epsDoppler = [obj.epsDoppler,eps_doppler];
            M=5; %Number of samples to average

            %{
            if(size(obj.epsDoppler,2)>M && eps_doppler > 1 && eps_doppler>mean(obj.epsDoppler(end-M:end-1)))
                %obj.epsDoppler =obj.epsDoppler(end-M:end-1);
                diffs = (obj.particles(:, [1 3])' - z)';
                likelihood_x = exp(-0.5 * (diffs(:, 1).^2) / obj.std_meas(1)^2);
                likelihood_y = exp(-0.5 * (diffs(:, 2).^2) / obj.std_meas(2)^(2));
                likelihood = likelihood_x .* likelihood_y;
                
                % Normalize the likelihood
                likelihood = likelihood / sum(likelihood);
                
                 weights_adapt =obj.weights .*likelihood;

                % Resample if too few effective particles
                neff = obj.NEFF(weights_adapt);
                if neff < obj.N / 2
                    indexes = obj.resampleSystematic(weights_adapt);
                    [obj.particles, weights_adapt] = obj.resampleFromIndex(obj.particles, indexes);
                end

                [meanValue,~] = obj.estimate(obj.particles, weights_adapt);
                
            else
            %}
                 
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
            %end

            X_est = meanValue;            
            PF_obj = obj;
        end
    end
    methods(Static)
        function [ resample_idx ] = resampleSystematic( w )
            N = length(w);
            % Stratified resampling
            resample_idx = zeros(N, 1);
            cumulative_weights = cumsum(w);
            thresholds = linspace(0, 1 - 1/N, N) + rand(1, N)/N;
        
            j = 1;
            for i = 1:N
                while thresholds(i) > cumulative_weights(j)
                    j = j + 1;
                end
                resample_idx(i) = j;
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