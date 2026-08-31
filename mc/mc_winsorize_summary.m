function mc_winsorize_summary(pct_lo, pct_hi)
% MC_WINSORIZE_SUMMARY  Recompute filter_perf mean+/-std with winsorization.
%
% Loads mc_results.mat and, for each (scenario, filter, metric) triple,
% clips the per-run values to the [pct_lo, pct_hi] percentile band across
% the N_MC runs before computing mean+/-std. This removes pathological-run
% dominance without dropping any run count from the reported statistic.
%
% Defaults: pct_lo=5, pct_hi=95.

    if nargin < 1 || isempty(pct_lo), pct_lo = 5;  end
    if nargin < 2 || isempty(pct_hi), pct_hi = 95; end

    S = load('mc_results.mat');
    raw = S.raw;                                  % N_MC x scen x filt x metric
    [Nmc, Nscen, Nfilt, Nmet] = size(raw);

    mu_w    = nan(Nscen, Nfilt, Nmet);
    sigma_w = nan(Nscen, Nfilt, Nmet);
    med_w   = nan(Nscen, Nfilt, Nmet);

    for si = 1:Nscen
        for fi = 1:Nfilt
            for mi = 1:Nmet
                v = squeeze(raw(:, si, fi, mi));
                v = v(~isnan(v));
                if isempty(v)
                    continue;
                end
                lo = prctile(v, pct_lo);
                hi = prctile(v, pct_hi);
                vw = min(max(v, lo), hi);
                mu_w(si, fi, mi)    = mean(vw);
                sigma_w(si, fi, mi) = std(vw);
                med_w(si, fi, mi)   = median(v);   % raw median for reference
            end
        end
    end

    save('mc_results_winsorized.mat', 'mu_w', 'sigma_w', 'med_w', ...
         'pct_lo', 'pct_hi');

    fprintf('\n=== Table filter_perf, winsorized at [%d, %d] pct ===\n', ...
            pct_lo, pct_hi);
    for si = 1:Nscen
        fprintf('\n[%s]\n', S.scenarios{si});
        fprintf('%-10s', 'Metric');
        for fi = 1:Nfilt, fprintf(' %-22s', S.filters{fi}); end
        fprintf('\n');
        for mi = 1:Nmet
            fprintf('%-10s', S.metrics{mi});
            for fi = 1:Nfilt
                fprintf(' %9.3f +/- %6.3f (med %6.3f)  ', ...
                        mu_w(si,fi,mi), sigma_w(si,fi,mi), med_w(si,fi,mi));
            end
            fprintf('\n');
        end
    end
end
