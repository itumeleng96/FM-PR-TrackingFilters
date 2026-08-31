% Test CSKF specifically with per-axis adaptive Q vs uniform Q vs no adaptive Q.
% FM PR emphasis: doppler is the fine-resolution axis; we care most about
% doppler RMSE.

clc; clear;

addpath('07_Evaluation/GroundTruth/', '06_Tracking/Filters/KF/');

scenarios = {'levelFlight', 'landing', 'takeoff', 'orbit360'};

dt = 1;
kf_std_meas = [4.9038, 0.9985];
kf_std_acc  = [0.0048354, 0.0991];

% FM PR resolution (for normalisation of the composite score):
d_range_km  = 1.5;   % ~c/fs/1000 = 1.5 km bin
d_dopp_hz   = 1.0;   % ~fs/N = 1 Hz per bin @ 1s integration

n_scen = numel(scenarios);
labels = {'Q_OFF', 'Q_ON_per_axis'};

fprintf('\n%-14s | ', 'Scenario');
for L = labels
    fprintf('%-22s | ', sprintf('%s (r/d)', L{1}));
end
fprintf('\n%s\n', repmat('-', 1, 90));

metrics = struct();
for s = 1:n_scen
    name = scenarios{s};
    cache = sprintf('meas_%s.mat', name);
    S = load(cache); d = S.data;

    fprintf('%-14s | ', name);
    for L = 1:numel(labels)
        adapt_q = strcmp(labels{L}, 'Q_ON_per_axis');
        initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
        flt = CSKF(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), initial_state);
        flt.cs_adapt_q = adapt_q;

        est_r = zeros(1, d.N); est_d = zeros(1, d.N);
        for k = 1:d.N
            [~,  flt] = flt.predict();
            [xk, flt] = flt.update(d.meas(:,k));
            est_r(k) = xk(1); est_d(k) = xk(3);
        end
        N = d.N_align;
        er = est_r(1:N)' - d.true_r;
        ed = est_d(1:N)' - d.true_d;
        rmse_r = sqrt(mean(er.^2));
        rmse_d = sqrt(mean(ed.^2));
        % FM-PR composite: bins-of-measurement-resolution equivalent
        score  = rmse_r/d_range_km + rmse_d/d_dopp_hz;

        fprintf('%6.2fkm/%6.1fHz s=%5.1f | ', rmse_r, rmse_d, score);
        metrics.(name).(labels{L}) = [rmse_r, rmse_d, score];
    end
    fprintf('\n');
end

fprintf('\nComposite score = RMSE_r/1.5km + RMSE_d/1Hz (bins-of-measurement units).\n');
fprintf('Doppler dominates because 1 Hz is the fine axis; range 1.5 km is the coarse axis.\n');

save('cskf_per_axis_q_test.mat', 'metrics');
