% =============================================================
% sanity_check_covscaling.m
%
% Purpose:
%   Confirm the refactored CSKF / CSUKF / CSRGNF still execute
%   after the eps^2 -> ek^2 swap in the R-blend line.
%   Injects synthetic measurements (constant-velocity target
%   + Gaussian noise + a couple of outliers) and runs a short
%   predict/update loop for each filter. Reports pass/fail and
%   a few state snapshots.
% =============================================================

clc; clear;

addpath('06_Tracking/Filters/KF/', ...
        '06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/RGNF/', ...
        '06_Tracking/Filters/PF/');

dt = 1;
N  = 40;
X0 = [1000; 5; 200; 0.2];               % [range, r_dot, doppler, dopp_dot]

% Per-filter parameters mirror runComputationalLoadRawFilter.m (cases 1/3/5/7)
kf_std_meas   = [4.9038, 0.9985];
kf_std_acc    = [0.0048354, 0.0991];
ukf_std_meas  = [0.9707, 0.79739];
ukf_std_acc   = [0.0076533, 0.09938];
ukf_alpha = 1e-4; ukf_kappa = 0; ukf_beta = 2;
pf_std_meas   = [10, 2];                % stdev (squared in likelihood)
pf_std_acc    = [1.429, 1.9452];
N_particles   = 2000;                   % lighter than 10000 for a sanity check
rgnf_std_meas = [2.046, 0.98];
rgnf_std_acc  = [0.057027, 0.047789];
rgnf_max_iter = 100; rgnf_lambda = 1;

% -------- synthetic measurements (constant-velocity + Gaussian noise) --------
F = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
truth = zeros(4, N);
truth(:,1) = X0;
for k = 2:N
    truth(:,k) = F*truth(:,k-1);
end
% Noise sized on the tightest R (KF) so all filters see the same measurement realisation.
z = [truth(1,:); truth(3,:)] + [sqrt(kf_std_meas(1)); sqrt(kf_std_meas(2))] .* randn(2, N);

% Inject two outliers to trigger the covariance-scaling branch
z(:,20) = z(:,20) + [30; 8];
z(:,30) = z(:,30) + [40; 10];

filters = { ...
    'CSKF',   @() CSKF(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), X0); ...
    'CSUKF',  @() CSUKF(dt, ukf_std_acc, ukf_std_meas(1), ukf_std_meas(2), X0', ukf_alpha, ukf_kappa, ukf_beta); ...
    'CSRGNF', @() CSRGNF(dt, rgnf_std_acc, rgnf_std_meas(1), rgnf_std_meas(2), X0, rgnf_max_iter, rgnf_lambda); ...
    'CSPF',   @() CSPF(dt, pf_std_acc, pf_std_meas, X0, N_particles); ...
};

fprintf('\n== covariance-scaling sanity check ==\n');
for f = 1:size(filters,1)
    name = filters{f,1};
    ctor = filters{f,2};
    try
        flt = ctor();
        est = zeros(4, N);
        for k = 1:N
            [~,   flt] = flt.predict();
            [xk,  flt] = flt.update(z(:,k));
            est(:,k) = xk(:);
        end
        err = est([1 3],:) - truth([1 3],:);
        rmse = sqrt(mean(err.^2, 2));
        fprintf('  %-7s  OK   RMSE range=%.2f  doppler=%.2f\n', name, rmse(1), rmse(2));
    catch ME
        fprintf('  %-7s  FAIL %s\n', name, ME.message);
    end
end
fprintf('=====================================\n');
