% Test adaptive Q by re-running per-scenario tuning across all four filters.
% Uses cached measurements. Reports RMSE with adaptive Q ON vs OFF for each
% (filter, scenario) pair.

clc; clear;

addpath('07_Evaluation/GroundTruth/', ...
        '06_Tracking/Filters/KF/', '06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/RGNF/', '06_Tracking/Filters/PF/');

scenarios = {'levelFlight', 'landing', 'takeoff', 'orbit360'};

dt = 1;
kf_std_meas   = [4.9038, 0.9985];   kf_std_acc   = [0.0048354, 0.0991];
ukf_std_meas  = [0.9707, 0.79739];  ukf_std_acc  = [0.0076533, 0.09938];
ukf_alpha=1e-4; ukf_kappa=0; ukf_beta=2;
pf_std_meas   = [10, 2];            pf_std_acc   = [1.429, 1.9452]; N_pf = 3000;
rgnf_std_meas = [2.046, 0.98];      rgnf_std_acc = [0.057027, 0.047789];
rgnf_max_iter = 100; rgnf_lambda = 1;

filters = {'CSKF','CSUKF','CSRGNF','CSPF'};
n_scen  = numel(scenarios);
n_filt  = numel(filters);

results = struct();

for adapt_q_flag = [false true]
    tag = char('QoffQon' * (adapt_q_flag+0) + 0);   % just for labelling
    if adapt_q_flag, tag = 'adaptQ_ON'; else, tag = 'adaptQ_OFF'; end
    fprintf('\n==== %s ====\n', tag);
    fprintf('%-14s | %-8s %-10s | %-8s %-10s | %-8s %-10s | %-8s %-10s\n', ...
            'Scenario','CSKF_r','CSKF_d','UKF_r','UKF_d','RGN_r','RGN_d','PF_r','PF_d');
    for s = 1:n_scen
        name = scenarios{s};
        cache = sprintf('meas_%s.mat', name);
        S = load(cache); d = S.data;
        row = zeros(1, n_filt*2);
        for f = 1:n_filt
            initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
            switch filters{f}
                case 'CSKF'
                    flt = CSKF(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), initial_state);
                case 'CSUKF'
                    flt = CSUKF(dt, ukf_std_acc, ukf_std_meas(1), ukf_std_meas(2), initial_state', ukf_alpha, ukf_kappa, ukf_beta);
                case 'CSRGNF'
                    flt = CSRGNF(dt, rgnf_std_acc, rgnf_std_meas(1), rgnf_std_meas(2), initial_state, rgnf_max_iter, rgnf_lambda);
                case 'CSPF'
                    flt = CSPF(dt, pf_std_acc, pf_std_meas, initial_state, N_pf);
            end
            if isprop(flt, 'cs_adapt_q'); flt.cs_adapt_q = adapt_q_flag; end

            est_r = zeros(1, d.N); est_d = zeros(1, d.N);
            for k = 1:d.N
                [~,  flt] = flt.predict();
                [xk, flt] = flt.update(d.meas(:,k));
                est_r(k) = xk(1); est_d(k) = xk(3);
            end
            N = d.N_align;
            er = est_r(1:N)' - d.true_r; ed = est_d(1:N)' - d.true_d;
            row(2*f-1) = sqrt(mean(er.^2));
            row(2*f)   = sqrt(mean(ed.^2));
        end
        fprintf('%-14s | %6.2fkm %6.1fHz | %6.2fkm %6.1fHz | %6.2fkm %6.1fHz | %6.2fkm %6.1fHz\n', ...
                name, row(1),row(2),row(3),row(4),row(5),row(6),row(7),row(8));
        results.(tag).(name) = row;
    end
end

save('adaptive_q_test.mat', 'results');
fprintf('\nSaved to adaptive_q_test.mat\n');
