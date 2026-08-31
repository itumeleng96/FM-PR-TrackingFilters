% Summarize scaling results for Table III.
% Computes per-scan total cost (update + predict) with proper Option B
% aggregation: for each MC run, sum update+predict per scan and take the
% per-run mean, then report mean +/- std across the 1000 per-run means.

clc; clear;

load('evaluationScripts/ComputationalLoad/computational_load_scaling_results.mat', 'results', 'warmup_scans');

N_labels = {'N=1', 'N=3', 'N=5', 'N=10'};
filter_labels = {'KF', 'UKF', 'PF', 'RGNF'};

fprintf('\n=======================================================================\n');
fprintf(' TABLE III (Option B): per-scan total cost mean +/- std (us)\n');
fprintf(' Aggregation: per-run mean of (update+predict) across scans, then\n');
fprintf(' mean +/- std across the 1000 per-run means.\n');
fprintf('=======================================================================\n');
fprintf('%-6s | %-15s | %-15s | %-15s | %-15s\n', 'Filter', N_labels{:});
fprintf('-----------------------------------------------------------------------\n');

for f = 1:numel(filter_labels)
    fname = filter_labels{f};
    row = sprintf('%-6s |', fname);
    for s = 1:numel(results)
        % results(s).r holds the mu/sd summary; the scaling script did NOT
        % save the raw timing matrices, so we fall back to an independent-sum
        % approximation using the reported update/predict mu/sd (correlation
        % unknown).
        %
        % Independent-sum approximation:
        %   mu_tot   = mu_update + mu_predict
        %   sigma_tot = sqrt(sigma_update^2 + sigma_predict^2)
        r = results(s).r;
        switch fname
            case 'KF'
                mu  = r.mu_uKF  + r.mu_pKF;
                sd  = sqrt(r.sd_uKF^2  + r.sd_pKF^2);
            case 'UKF'
                mu  = r.mu_uUKF + r.mu_pUKF;
                sd  = sqrt(r.sd_uUKF^2 + r.sd_pUKF^2);
            case 'PF'
                mu  = r.mu_uPF  + r.mu_pPF;
                sd  = sqrt(r.sd_uPF^2  + r.sd_pPF^2);
            case 'RGNF'
                mu  = r.mu_uRG  + r.mu_pRG;
                sd  = sqrt(r.sd_uRG^2  + r.sd_pRG^2);
        end
        row = [row, sprintf(' %7.0f +/- %5.0f |', mu, sd)]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\nAlso reporting scaling ratios (N=10 total / N=1 total):\n');
for f = 1:numel(filter_labels)
    fname = filter_labels{f};
    r1  = results(1).r;
    r10 = results(4).r;
    switch fname
        case 'KF';   tot1 = r1.mu_uKF  + r1.mu_pKF;   tot10 = r10.mu_uKF  + r10.mu_pKF;
        case 'UKF';  tot1 = r1.mu_uUKF + r1.mu_pUKF;  tot10 = r10.mu_uUKF + r10.mu_pUKF;
        case 'PF';   tot1 = r1.mu_uPF  + r1.mu_pPF;   tot10 = r10.mu_uPF  + r10.mu_pPF;
        case 'RGNF'; tot1 = r1.mu_uRG  + r1.mu_pRG;   tot10 = r10.mu_uRG  + r10.mu_pRG;
    end
    fprintf('  %-4s: cost_ratio = %.2fx\n', fname, tot10/tot1);
end
