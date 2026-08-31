% Compare cached measurements against ground truth for each scenario.
% If measurement residual std is comparable to or larger than filter RMSE,
% the pipeline (CFAR/mean-shift) is the bottleneck — no filter tuning
% will help.

clc; close all;

addpath('07_Evaluation/GroundTruth/');
scenarios = {'levelFlight', 'landing', 'takeoff', 'orbit360'};

fprintf('\n%-14s | %8s %8s | %8s %8s | %8s %8s\n', ...
        'Scenario','res_r_mn','res_r_sd','res_d_mn','res_d_sd','miss%_r','miss%_d');
fprintf('%s\n', repmat('-', 1, 80));

for s = 1:numel(scenarios)
    name = scenarios{s};
    cache = sprintf('meas_%s.mat', name);
    S = load(cache); d = S.data;
    N = d.N_align;

    meas_r = d.meas(1, 1:N)';
    meas_d = d.meas(2, 1:N)';

    res_r = meas_r - d.true_r;
    res_d = meas_d - d.true_d;

    % "Missed" detections: measurement within a bin of zero (no CFAR
    % detection, filled with prior scan or zero).
    miss_r = mean(abs(meas_r) < 0.1) * 100;
    miss_d = mean(abs(meas_d) < 0.1) * 100;

    fprintf('%-14s | %8.2f %8.2f | %8.2f %8.2f | %7.1f%% %7.1f%%\n', ...
            name, mean(res_r), std(res_r), mean(res_d), std(res_d), ...
            miss_r, miss_d);
end
