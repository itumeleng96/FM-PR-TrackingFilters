function [t, bistatic_range_km, doppler_hz, range_rate_kmps, doppler_rate_hzps] = computeGroundTruth(scenario)
% COMPUTEGROUNDTRUTH  Bistatic range/Doppler ground truth per scenario.
%
%   [t, R_km, fd_Hz] = computeGroundTruth('levelFlight')
%
% Scenarios (matching the paper):
%   'levelFlight'  - Single-target level flight (60 s, 4 waypoints)
%   'landing'      - Landing manoeuvre       (120 s, 5 waypoints)
%   'takeoff'      - Take-off manoeuvre      ( 60 s, 9 waypoints)
%   'orbit360'     - 360-degree orbit        (120 s, 6 waypoints)
%
% Waypoints are interpolated with a variational (natural) cubic spline
% (csape) and sampled at 1-s intervals; range and Doppler are then
% computed from the bistatic geometry using the shared Tx/Rx placement.
%
% Bistatic formulas (Griffiths & Baker 2005, Bistatic Radar Systems, Ch 3):
%   Bistatic range   R_b = |p-Tx| + |p-RefRx| - Baseline
%   Bistatic Doppler f_d = -(1/lambda) d/dt (|p-Tx| + |p-RefRx|)
%
% Sign convention: approaching target -> f_d > 0 (matches
% calculateBistaticCubicFinal.m and existing evaluation scripts).

    % ---- Bistatic geometry (shared across scenarios) ----
    Tx_Pos     = [ 6440;  10760;  1056];
    RefRx_Pos  = [    0;      0;  1000];
    Baseline_m = norm(RefRx_Pos - Tx_Pos);
    fc         = 94e6;                          % FM carrier
    c          = 299792458;
    lambda     = c / fc;

    % ---- Per-scenario waypoints (metres, seconds) ----
    switch lower(scenario)
        case 'levelflight'
            TargetPos   = [ 15000, 22000, 3600;
                             3000, 21000, 3600;
                            -2000, 21000, 3600;
                            -6000, 22000, 3600];
            TargetTimes = [0, 30, 45, 60];
        case 'landing'
            TargetPos   = [ -5400,  10000, 5000;
                           -15000,   6000, 3000;
                           -18000,      0, 2000;
                           -15000,  -6000, 1000;
                            -8500,  -5000,    0];
            TargetTimes = [0, 30, 60, 90, 120];
        case 'takeoff'
            TargetPos   = [ -7800,  -2000, 1000;
                            -8200,      0, 1800;
                            -8500,   1500, 2200;
                            -9000,   3000, 2800;
                            -9500,   4500, 3500;
                           -10000,   6000, 4200;
                           -10500,   7500, 4800;
                           -11000,   9000, 5400;
                           -11500,  10500, 6000];
            TargetTimes = [0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60];
        case 'orbit360'
            TargetPos   = [  8000, 31000, 8000;
                              -528, 30851, 8000;
                             -4354, 34339, 8000;
                                 3, 34783, 8000;
                             -5311, 29938, 8000;
                            -11431, 27706, 8000];
            TargetTimes = [0, 35, 50, 75, 95, 120];
        case '3targets'
            % Three straight-line targets from FERS/BackupScenarios/scenario_N3_targets.fersxml.
            % Each target is one 2-waypoint linear path over 0-60 s.
            multi_targets = struct( ...
                'T1', struct('pos', [13000, 12000, 3000; 13000, 18000, 3000], 'times', [0, 60]), ...
                'T5', struct('pos', [22000, 15000, 3600; 22000, 21000, 3600], 'times', [0, 60]), ...
                'T10', struct('pos', [35000, 12000, 3200; 35000,  6000, 3200], 'times', [0, 60]));
            t = (0:1:60)';
            names = fieldnames(multi_targets);
            NT = numel(names);
            bistatic_range_km = zeros(numel(t), NT);
            doppler_hz        = zeros(numel(t), NT);
            for ti = 1:NT
                mt = multi_targets.(names{ti});
                pos_t = interp1(mt.times, mt.pos, t, 'linear').';   % 3 x N (linear for 2-point paths)
                r_tx  = vecnorm(pos_t - Tx_Pos);
                r_ref = vecnorm(pos_t - RefRx_Pos);
                total_path_m = r_tx + r_ref;
                bistatic_range_km(:, ti) = (total_path_m - Baseline_m).' / 1000;
                d_path = zeros(1, numel(t));
                d_path(1) = total_path_m(2) - total_path_m(1);
                d_path(end) = total_path_m(end) - total_path_m(end-1);
                if numel(t) > 2
                    d_path(2:end-1) = (total_path_m(3:end) - total_path_m(1:end-2)) / 2;
                end
                doppler_hz(:, ti) = (-d_path ./ (1 * lambda)).';
            end
            return;
        otherwise
            error('computeGroundTruth:unknownScenario', ...
                  'Unknown scenario "%s". Use levelFlight | landing | takeoff | orbit360 | 3targets.', scenario);
    end

    % ---- Cubic-spline interpolation at 1 s resolution ----
    t         = (TargetTimes(1) : 1 : TargetTimes(end))';
    spline_x  = csape(TargetTimes, TargetPos(:,1), 'variational');
    spline_y  = csape(TargetTimes, TargetPos(:,2), 'variational');
    spline_z  = csape(TargetTimes, TargetPos(:,3), 'variational');
    pos       = [fnval(spline_x, t)'; fnval(spline_y, t)'; fnval(spline_z, t)'];   % 3xN

    % ---- Bistatic range (km) ----
    r_tx  = vecnorm(pos - Tx_Pos);
    r_ref = vecnorm(pos - RefRx_Pos);
    total_path_m       = r_tx + r_ref;                     % 1xN (m)
    bistatic_range_km  = (total_path_m - Baseline_m).' / 1000;

    % ---- Bistatic Doppler (Hz) via central difference on total path ----
    % f_d = -(1/lambda) d(total_path)/dt   (positive when path shortens)
    dt = 1;             % s
    N  = numel(t);
    d_path = zeros(1, N);
    d_path(1)     = total_path_m(2) - total_path_m(1);
    d_path(N)     = total_path_m(N) - total_path_m(N-1);
    if N > 2
        d_path(2:N-1) = (total_path_m(3:N) - total_path_m(1:N-2)) / 2;
    end
    doppler_hz = (-d_path ./ (dt * lambda)).';
end
