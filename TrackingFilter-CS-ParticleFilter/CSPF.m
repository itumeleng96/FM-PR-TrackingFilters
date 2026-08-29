classdef CSPF

    properties
        dt,             %%Sampling Time
        F,              %The state transition matrix
        H,              %Measurement matrix (2x4, maps state to [range; doppler])
        X,              %Combined state estimate (weighted particle mean, 4x1)
        Q,              %The Process Noise Covariance Matrix
        particles,      %Matrix containing the current State of the Particles
        weights,        %A vector containig the current weights of the particles
        N,              %Number of particles to use 
        scaling_factor; % 
        std_meas;
        R;              % Adaptive measurement noise covariance (2x2 diagonal)
        R_nominal;      % Nominal (initial) R
        Q_nominal;      % Nominal (initial) Q
        S;
        P;
        cs_alpha;       % Akhlaghi forgetting factor for R
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

            obj.H = [1 0 0 0;
                     0 0 1 0];
                    
           
            

            eps = 1e-7;

            obj.Q = [std_acc(1)*(dt^4)/4 + eps, std_acc(1)*(dt^3)/2, 0, 0;
                     std_acc(1)*(dt^3)/2, std_acc(1)*dt^2 + eps, 0, 0;
                     0, 0, std_acc(2)*(dt^4)/4 + eps, std_acc(2)*(dt^3)/2;
                     0, 0, std_acc(2)*(dt^3)/2, std_acc(2)*dt^2 + eps];

            
            % Measurement noise covariance (persistent across scans; adapted
            % per Eq. (9a) on outlier trigger, relaxed per Eq. (9b) otherwise).
            obj.R         = diag(std_meas(:).^2);
            obj.R_nominal = obj.R;
            obj.Q_nominal = obj.Q;

            % Akhlaghi (2018) forgetting factor for R adaptation.
            obj.cs_alpha  = 0.98;
        end
        
        function [X_pred, PF_obj] = predict(obj)
            % SIR bootstrap prediction: draw x_k^i ~ p(x_k | x_{k-1}^i).
            %
            % Discrete white-noise-acceleration transition model:
            %   x_k^i = F x_{k-1}^i + w_k^i,   w_k^i ~ N(0, Q)
            %
            % Predicted state estimate is the weighted empirical mean.
            % Innovation covariance follows the standard Kalman form:
            %   S_k = H P_k^- H^T + R
            %
            % Refs: Gordon, Salmond & Smith (1993), "Novel approach to
            %       nonlinear/non-Gaussian Bayesian state estimation";
            %       Bar-Shalom, Li & Kirubarajan (2001), "Estimation with
            %       Applications to Tracking and Navigation", Sec. 5.2.

            obj.particles = obj.particles * obj.F' + randn(obj.N, 4) * chol(obj.Q);

            [meanVal, varVal] = obj.estimate(obj.particles, obj.weights);
            X_pred = meanVal(:);
            obj.X  = X_pred;

            obj.P = diag(varVal);
            obj.S = diag([obj.P(1,1); obj.P(3,3)]) + obj.R;

            PF_obj = obj;
        end
        
        function [X_est, PF_obj] = update(obj, z)
            % Akhlaghi (2018) continuous R update in PF form. Q is used at
            % prediction but not adapted online: Akhlaghi's Q update presumes
            % an explicit Kalman gain, which the PF does not have.
            %   R <- alpha R + (1-alpha)( eps eps' + H P^- H' )
            %   H P^- H' is the empirical measurement covariance of the
            %   predicted particle cloud.
            % Outlier rejection is done upstream at the association step.

            z = z(:);
            alphaFactor = obj.cs_alpha;

            % Pre-update state and predicted particle-cloud measurement covariance.
            HPH = diag([obj.P(1,1); obj.P(3,3)]);

            % SIR weight update using current R for the likelihood.
            sigma2_x = obj.R(1,1);
            sigma2_y = obj.R(2,2);
            diffs = obj.particles(:, [1 3]) - z.';
            log_L = -0.5 * ( diffs(:,1).^2 ./ sigma2_x + diffs(:,2).^2 ./ sigma2_y );
            log_L = log_L - max(log_L);
            L     = exp(log_L);

            obj.weights = obj.weights .* L + realmin;
            obj.weights = obj.weights ./ sum(obj.weights);

            if obj.NEFF(obj.weights) < obj.N / 2
                indexes = obj.resampleSystematic(obj.weights);
                [obj.particles, obj.weights] = obj.resampleFromIndex(obj.particles, indexes);
            end

            [meanValPost, ~] = obj.estimate(obj.particles, obj.weights);
            X_est   = meanValPost(:);
            obj.X   = X_est;

            % Akhlaghi Eq. 11 — continuous R update, PF form.
            eps_post = z - meanValPost([1; 3]);
            obj.R    = alphaFactor * obj.R + (1 - alphaFactor) * ...
                       (eps_post * eps_post.' + HPH);

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
        

        function [meanVal, varVal] = estimate(particles, weights)
            % Weighted empirical moments of the particle cloud.
            %
            %   E[x]   = sum_i w_i x_i / sum_i w_i
            %   Var[x] = sum_i w_i (x_i - E[x])^2 / sum_i w_i
            %
            % Ref: Ristic, Arulampalam & Gordon (2004),
            %      "Beyond the Kalman Filter", Sec. 3.3.
            w       = weights(:);
            wsum    = sum(w);
            meanVal = (particles' * w) ./ wsum;             % 4x1
            dev     = particles - meanVal.';                % Nx4
            varVal  = (dev.^2)' * w ./ wsum;                % 4x1
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