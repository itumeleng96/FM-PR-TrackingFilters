classdef IMM_CV2
    % Two-model IMM: CV_low (cruise) and CV_high (manoeuvre), same state
    % dim so no augmented-state mixing needed.
    % Refs: Bar-Shalom, Li & Kirubarajan (2001), "Estimation with Applications
    %       to Tracking and Navigation", Ch. 11.6.

    properties
        dt, F, H, R
        Q_low, Q_high            % per-model process-noise covariances
        X1, P1, X2, P2           % per-model state and covariance
        mu                        % 2x1 model probabilities
        TPM                       % 2x2 Markov transition matrix
        X, P, S                   % combined output (compat with runner)
        residuals_x, residuals_y  % legacy fields for runner trigger detection
    end

    methods
        function obj = IMM_CV2(dt, std_acc_low, std_acc_high, r_std, rdot_std, X_initial)
            obj.dt = dt;
            obj.F = [1 dt 0 0;
                     0  1 0 0;
                     0  0 1 dt;
                     0  0 0 1];
            obj.H = [1 0 0 0;
                     0 0 1 0];
            obj.R = [r_std 0; 0 rdot_std];

            % Model 1: cruise (low process noise)
            obj.Q_low  = [std_acc_low(1)*dt^4/4, std_acc_low(1)*dt^3/2, 0, 0;
                          std_acc_low(1)*dt^3/2, std_acc_low(1)*dt^2,   0, 0;
                          0, 0, std_acc_low(2)*dt^4/4, std_acc_low(2)*dt^3/2;
                          0, 0, std_acc_low(2)*dt^3/2, std_acc_low(2)*dt^2];

            % Model 2: manoeuvre (high process noise)
            obj.Q_high = [std_acc_high(1)*dt^4/4, std_acc_high(1)*dt^3/2, 0, 0;
                          std_acc_high(1)*dt^3/2, std_acc_high(1)*dt^2,   0, 0;
                          0, 0, std_acc_high(2)*dt^4/4, std_acc_high(2)*dt^3/2;
                          0, 0, std_acc_high(2)*dt^3/2, std_acc_high(2)*dt^2];

            % Init: both models start at same state
            obj.X1 = X_initial;   obj.P1 = diag([5 1 2 1]);
            obj.X2 = X_initial;   obj.P2 = diag([5 1 2 1]);

            % Init model probs (equal), Markov transition (sticky)
            obj.mu  = [0.5; 0.5];
            obj.TPM = [0.95 0.05;
                       0.05 0.95];

            obj.X = X_initial;  obj.P = diag([5 1 2 1]);  obj.S = obj.R;

            obj.residuals_x = [];  obj.residuals_y = [];
        end

        function [X_pred, obj] = predict(obj)
            % IMM step 1: mixing probabilities
            c   = obj.TPM.' * obj.mu;          % 2x1  normalisation
            mix = zeros(2,2);                   % mix(i,j) = P(model i at k-1 | model j at k)
            for j = 1:2
                for i = 1:2
                    mix(i,j) = obj.TPM(i,j) * obj.mu(i) / max(c(j), 1e-12);
                end
            end

            % IMM step 2: mixed initial conditions per model
            X01 = mix(1,1)*obj.X1 + mix(2,1)*obj.X2;
            X02 = mix(1,2)*obj.X1 + mix(2,2)*obj.X2;
            d1a = obj.X1 - X01; d2a = obj.X2 - X01;
            d1b = obj.X1 - X02; d2b = obj.X2 - X02;
            P01 = mix(1,1)*(obj.P1 + d1a*d1a.') + mix(2,1)*(obj.P2 + d2a*d2a.');
            P02 = mix(1,2)*(obj.P1 + d1b*d1b.') + mix(2,2)*(obj.P2 + d2b*d2b.');

            % IMM step 3: model-conditioned prediction
            obj.X1 = obj.F * X01;   obj.P1 = obj.F * P01 * obj.F.' + obj.Q_low;
            obj.X2 = obj.F * X02;   obj.P2 = obj.F * P02 * obj.F.' + obj.Q_high;

            % Combined output (probability-weighted)
            obj.X = obj.mu(1)*obj.X1 + obj.mu(2)*obj.X2;
            d1 = obj.X1 - obj.X;  d2 = obj.X2 - obj.X;
            obj.P = obj.mu(1)*(obj.P1 + d1*d1.') + obj.mu(2)*(obj.P2 + d2*d2.');
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            X_pred = obj.X;
        end

        function [X_est, obj] = update(obj, z)
            % Per-model innovation, likelihood, KF update
            d1 = z - obj.H * obj.X1;   S1 = obj.H * obj.P1 * obj.H.' + obj.R;
            d2 = z - obj.H * obj.X2;   S2 = obj.H * obj.P2 * obj.H.' + obj.R;
            L1 = (2*pi)^(-1) / sqrt(max(det(S1),1e-12)) * exp(-0.5 * d1.' * (S1 \ d1));
            L2 = (2*pi)^(-1) / sqrt(max(det(S2),1e-12)) * exp(-0.5 * d2.' * (S2 \ d2));

            K1 = (obj.P1 * obj.H.') / S1;
            obj.X1 = obj.X1 + K1 * d1;
            obj.P1 = (eye(4) - K1 * obj.H) * obj.P1;

            K2 = (obj.P2 * obj.H.') / S2;
            obj.X2 = obj.X2 + K2 * d2;
            obj.P2 = (eye(4) - K2 * obj.H) * obj.P2;

            % IMM step 4: model probability update
            c    = obj.TPM.' * obj.mu;
            L    = [L1; L2];
            post = L .* c;
            post = post / max(sum(post), 1e-12);
            obj.mu = post;

            % Combined output
            obj.X = obj.mu(1)*obj.X1 + obj.mu(2)*obj.X2;
            e1 = obj.X1 - obj.X;  e2 = obj.X2 - obj.X;
            obj.P = obj.mu(1)*(obj.P1 + e1*e1.') + obj.mu(2)*(obj.P2 + e2*e2.');
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            X_est = obj.X;
        end
    end
end
