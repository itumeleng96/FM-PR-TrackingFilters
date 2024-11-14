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
        epsDoppler;
        epsRange;
        P;

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
                    
           
            

            eps = 1e-7;

            obj.Q = [std_acc(1)*(dt^4)/4 + eps, std_acc(1)*(dt^3)/2, 0, 0;
                     std_acc(1)*(dt^3)/2, std_acc(1)*dt^2 + eps, 0, 0;
                     0, 0, std_acc(2)*(dt^4)/4 + eps, std_acc(2)*(dt^3)/2;
                     0, 0, std_acc(2)*(dt^3)/2, std_acc(2)*dt^2 + eps];

            
            obj.count =0;
            obj.updater =0;
            obj.update1 =0;

            obj.epsDoppler =[];
            obj.epsRange =[];




        end
        
        function [X_pred, PF_obj] = predict(obj)

            % Generate random Gaussian noise with zero mean and covariance matrix Q
            noise = randn(obj.N, 4) * chol(obj.Q);


            % Add process noise to particle states
            obj.particles(:, 1:4) = (obj.F(1:4, 1:4) * obj.particles(:, 1:4)' + noise(:, 1:4)')';
                    
            X_pred = mean(obj.particles, 1)';
            

            [meanValueS,~] = obj.estimate(obj.particles, obj.weights);
            deviations = obj.particles(:, [1 3]) - meanValueS;  % Deviation of particles from mean (Nx2)
            weighted_deviations = deviations .* sqrt(obj.weights);  % Apply weights (element-wise multiplication)
            covariance_matrix = (weighted_deviations' * weighted_deviations) / sum(obj.weights);  % Weighted covariance
            obj.P = [covariance_matrix(1,1),0,0,0;
                     0,0,0,0;
                     0,0,covariance_matrix(2,2),0;
                     0,0,0,0];

            slikelihood= covariance_matrix+ obj.std_meas;
            obj.S = [mean(slikelihood(:,1)),0;0,mean(slikelihood(:,2));];
            PF_obj = obj;

        end
        
        function [X_est, PF_obj] = update(obj, z)
        
            [meanValueS, ~] = obj.estimate(obj.particles, obj.weights);
            deviations = obj.particles(:, [1 3]) - meanValueS;  % Deviations of particles from the mean (Nx2)
            weighted_deviations = deviations .* sqrt(obj.weights);  % Apply weights element-wise
            covariance_matrix = (weighted_deviations' * weighted_deviations) / sum(obj.weights);  % Weighted covariance
            
            % Update the covariance matrix (P)
            obj.P = [covariance_matrix(1,1), 0, 0, 0;
                     0, 0, 0, 0;
                     0, 0, covariance_matrix(2,2), 0;
                     0, 0, 0, 0];
        
            % Calculate S (measurement noise + covariance)
            slikelihood = covariance_matrix + obj.std_meas;
            obj.S = [mean(slikelihood(:, 1)), 0;
                     0, mean(slikelihood(:, 2))];
        
            % Calculate residual (measurement innovation)
            ek = z - meanValueS';
            eps_range = ek(1)^2 / obj.S(1,1);  % Mahalanobis distance for Range
            eps_doppler = ek(2)^2 / obj.S(2,2);  % Mahalanobis distance for Doppler
        
            % Moving windows for residuals (Range and Doppler)
            obj.epsRange = [obj.epsRange, eps_range];
            obj.epsDoppler = [obj.epsDoppler, eps_doppler];
            M = 6;  % Number of samples to average
            alphaFactor =0.8;
            % Flags to track if outliers were detected
            outlier_detected_range = false;
            outlier_detected_doppler = false;
        
            % Outlier rejection for Doppler
            if size(obj.epsDoppler, 2) > M && eps_doppler > 1 && eps_doppler > (std(obj.epsDoppler(end-M:end-1))+mean(obj.epsDoppler(end-M:end-1)))
                outlier_detected_doppler = true;
            end
        
            % Outlier rejection for Range (x component)
            if size(obj.epsRange, 2) > M && eps_range > 2 && eps_range > (std(obj.epsRange(end-M:end-1)) +mean(obj.epsRange(end-M:end-1)))
                outlier_detected_range = true;
            end
        
            % Update likelihood calculations based on outlier detection
            diffs = (obj.particles(:, [1 3])' - z)';
            
            % Handle outliers in the x (Range) axis
            if outlier_detected_range
                % If Range outlier detected, soften the likelihood by adjusting the variance using alphaFactor
                likelihood_x = exp(-0.5 * (diffs(:, 1).^2) / (alphaFactor * obj.std_meas(1)^2 + (1 - alphaFactor) * obj.P(1,1)));
            else
                % Regular likelihood calculation for x (Range)
                likelihood_x = exp(-0.5 * (diffs(:, 1).^2) / obj.std_meas(1)^2);
            end
        
            % Handle outliers in the y (Doppler) axis
            if outlier_detected_doppler
                % If Doppler outlier detected, soften the likelihood by adjusting the variance using alphaFactor
                likelihood_y = exp(-0.5 * (diffs(:, 2).^2) / (alphaFactor * obj.std_meas(2)^2 + (1 - alphaFactor) * obj.P(3,3)));
            else
                % Regular likelihood calculation for y (Doppler)
                likelihood_y = exp(-0.5 * (diffs(:, 2).^2) / obj.std_meas(2)^2);
            end
        
            % Combine the likelihoods
            likelihood = likelihood_x .* likelihood_y;
        
            % Normalize the likelihood
            likelihood = likelihood / sum(likelihood);
        
            % Update the particle weights
            obj.weights = obj.weights .* likelihood;
            obj.weights = obj.weights + 1.e-300;  % Avoid zero weights
            obj.weights = obj.weights / sum(obj.weights);
        
            % Resample if too few effective particles
            neff = obj.NEFF(obj.weights);
            if neff < obj.N / 2
                indexes = obj.resampleSystematic(obj.weights);
                [obj.particles, obj.weights] = obj.resampleFromIndex(obj.particles, indexes);
            end
        
            % Recalculate the state estimate
            [meanValue, ~] = obj.estimate(obj.particles, obj.weights);
        
            % Output the estimated state
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