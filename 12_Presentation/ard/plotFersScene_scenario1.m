% plotFersScene_scenario1.m
% ------------------------------------------------------------------
% Renders the FERS "scenario_1_singleFile" scene in 3-D for Slide 4 of
% the CIE 2026 presentation. Geometry and waypoints taken verbatim from:
%   MscEng/FM-PR-TrackingFilters/01_FERS/BackupScenarios/scenario_1_singleFile.fersxml
%
% Style matches calculateBistaticCubicFinal.m (plot3 + text labels + view(3)).
%
% Output:
%   Figures/fers_scene_scenario1.svg
%   Figures/fers_scene_scenario1.png   (300 dpi fallback)
%
% Author: Itumeleng Malemela
% ------------------------------------------------------------------

clc; clear; close all;

% -------------------- Geometry (km) ---------------------------
% Stations
Tx_Pos     = [6.440;  10.760;  1.056];   % PrieskaTx (94 MHz, isotropic)
RefRx_Pos  = [0;      0;       1.000];   % Reference Rx (Yagi, az 42 deg)
SurvRx_Pos = [0;      0.001;   1.000];   % Surveillance Rx (Yagi, az 180 deg)

Baseline = norm(RefRx_Pos - Tx_Pos);   % ~ 12.54 km

% Target waypoints from scenario_1_singleFile.fersxml (cubic interp)
TargetWayPoints = [0, 30, 45, 60];
TargetPos = [ 15,  22,  3.6;
               3,  21,  3.6;
              -2,  21,  3.6;
              -6,  22,  3.6];

% -------------------- Cubic-spline trajectory ---------------------
t_interp = linspace(TargetWayPoints(1), TargetWayPoints(end), 241);
spline_x = csape(TargetWayPoints, TargetPos(:,1), 'variational');
spline_y = csape(TargetWayPoints, TargetPos(:,2), 'variational');
spline_z = csape(TargetWayPoints, TargetPos(:,3), 'variational');
x_i = fnval(spline_x, t_interp);
y_i = fnval(spline_y, t_interp);
z_i = fnval(spline_z, t_interp);

% -------------------- Figure --------------------------------------
fig = figure('Color','w','Units','pixels','Position',[100 100 960 640]);
hold on; grid on; box on;

% ---- Modern colour palette ---------------------------------------
col.trajectory = [0.12 0.35 0.74];   % deep cobalt blue
col.baseline   = [0.35 0.35 0.42];   % soft charcoal
col.txPath     = [0.86 0.30 0.28];   % coral red
col.rxPath     = [0.18 0.55 0.78];   % steel blue
col.tx         = [0.86 0.30 0.28];
col.refrx      = [0.18 0.62 0.35];   % emerald green
col.survrx     = [0.20 0.44 0.78];
col.ground     = [0.94 0.95 0.97];   % near-white bluish tint

% ---- Ground plane (soft tinted panel + subtle grid) --------------
[gx, gy] = meshgrid(linspace(-8, 18, 2), linspace(-2, 24, 2));
gz = zeros(size(gx));
surf(gx, gy, gz, 'FaceColor', col.ground, 'EdgeColor','none', ...
     'FaceAlpha', 0.55);

% Baseline (Tx <-> RefRx) — dash-dot, medium weight
plot3([Tx_Pos(1) RefRx_Pos(1)], [Tx_Pos(2) RefRx_Pos(2)], ...
      [Tx_Pos(3) RefRx_Pos(3)], '-.', 'Color', col.baseline, 'LineWidth',1.6);

% Signal-path illustration at trajectory midpoint
mid_idx = round(numel(t_interp)/2);
target_mid = [x_i(mid_idx); y_i(mid_idx); z_i(mid_idx)];
plot3([Tx_Pos(1) target_mid(1)], [Tx_Pos(2) target_mid(2)], ...
      [Tx_Pos(3) target_mid(3)], '--', 'Color', col.txPath, 'LineWidth',1.6);
plot3([target_mid(1) SurvRx_Pos(1)], [target_mid(2) SurvRx_Pos(2)], ...
      [target_mid(3) SurvRx_Pos(3)], '--', 'Color', col.rxPath, 'LineWidth',1.6);

% ---- Target trajectory (single clean line) -----------------------
plot3(x_i, y_i, z_i, '-', 'Color', col.trajectory, 'LineWidth', 3.2);
% Waypoint markers
plot3(TargetPos(:,1), TargetPos(:,2), TargetPos(:,3), 'o', ...
      'MarkerFaceColor', col.trajectory, 'MarkerEdgeColor','w', ...
      'MarkerSize',6, 'LineWidth',1);
% Direction arrowhead at end of trajectory
plot3(x_i(end), y_i(end), z_i(end), '^', ...
      'MarkerFaceColor', col.trajectory, 'MarkerEdgeColor','w', ...
      'MarkerSize',10, 'LineWidth',1);

% -------------------- Airplane STL at start (t = 0 s) -------------
stlFilePath = fullfile(fileparts(mfilename('fullpath')), 'plane3d.stl');
if exist(stlFilePath, 'file')
    airplane_model = stlread(stlFilePath);
    V_pl = airplane_model.Points;
    F_pl = airplane_model.ConnectivityList;

    % --- Dynamic scaling: fit the plane to a target length in scene ---
    bbox    = max(V_pl) - min(V_pl);
    nativeL = max(bbox);                    % longest STL dimension
    targetL = 3;                            % desired plane length in km
    V_pl    = V_pl * (targetL / nativeL);

    % Recentre STL at its own origin before rotating
    V_pl = V_pl - mean(V_pl, 1);

    % Yaw to align nose along the local trajectory heading at the
    % midpoint (where the Tx- and Rx-side signal paths converge).
    % NOTE: plane3d.stl's native nose axis is not +X. Tweak stlYawOffset
    % if the aircraft points sideways / backwards:
    %   stlYawOffset =   0     if native nose is along +X
    %   stlYawOffset = -pi/2   if native nose is along +Y
    %   stlYawOffset =  pi/2   if native nose is along -Y
    %   stlYawOffset =  pi     if native nose is along -X
    stlYawOffset = pi/2;
    dx = x_i(mid_idx+1) - x_i(mid_idx-1);
    dy = y_i(mid_idx+1) - y_i(mid_idx-1);
    heading = atan2(dy, dx) + stlYawOffset;
    R = [cos(heading), -sin(heading), 0;
         sin(heading),  cos(heading), 0;
         0,             0,            1];
    V_pl = (R * V_pl')';

    % Translate to the signal-path convergence point (trajectory midpoint)
    V_pl = V_pl + target_mid(:)';

    patch('Vertices', V_pl, 'Faces', F_pl, ...
          'FaceColor', [0.80 0.80 1.00], 'EdgeColor','none', ...
          'FaceLighting','gouraud');
    camlight('headlight');
    lighting gouraud;
else
    plot3(x_i(1), y_i(1), z_i(1), 'ko', ...
          'MarkerFaceColor','k', 'MarkerSize',9);
    warning('plane3d.stl not found at %s — using circle marker instead.', stlFilePath);
end

% -------------------- Station markers -----------------------------
plot3(Tx_Pos(1), Tx_Pos(2), Tx_Pos(3), '^', ...
      'MarkerFaceColor', col.tx, 'MarkerEdgeColor','w', ...
      'MarkerSize',13, 'LineWidth',1.2);
plot3(RefRx_Pos(1), RefRx_Pos(2), RefRx_Pos(3), 's', ...
      'MarkerFaceColor', col.refrx, 'MarkerEdgeColor','w', ...
      'MarkerSize',13, 'LineWidth',1.2);
plot3(SurvRx_Pos(1), SurvRx_Pos(2), SurvRx_Pos(3), 's', ...
      'MarkerFaceColor', col.survrx, 'MarkerEdgeColor','w', ...
      'MarkerSize',13, 'LineWidth',1.2);

% Ground drop-lines
plot3([Tx_Pos(1) Tx_Pos(1)],       [Tx_Pos(2) Tx_Pos(2)],       [0 Tx_Pos(3)],    ':', 'Color', col.tx,    'LineWidth',1.0);
plot3([RefRx_Pos(1) RefRx_Pos(1)], [RefRx_Pos(2) RefRx_Pos(2)], [0 RefRx_Pos(3)], ':', 'Color', col.refrx, 'LineWidth',1.0);
plot3([SurvRx_Pos(1) SurvRx_Pos(1)],[SurvRx_Pos(2) SurvRx_Pos(2)],[0 SurvRx_Pos(3)],':','Color', col.survrx,'LineWidth',1.0);

% -------------------- Labels --------------------------------------
labelFont = 'Helvetica';
text(Tx_Pos(1)+0.4, Tx_Pos(2)+0.4, Tx_Pos(3)+0.5, ...
     'Tx  (94 MHz)', 'FontSize',11, 'FontWeight','bold', ...
     'FontName',labelFont, 'Color', col.tx);
text(RefRx_Pos(1)-3.2, RefRx_Pos(2)-1.5, RefRx_Pos(3)+0.4, ...
     'Reference Rx', 'FontSize',11, 'FontWeight','bold', ...
     'FontName',labelFont, 'Color', col.refrx);
text(SurvRx_Pos(1)-3.5, SurvRx_Pos(2)-1.5, SurvRx_Pos(3)-0.4, ...
     'Surveillance Rx', 'FontSize',11, 'FontWeight','bold', ...
     'FontName',labelFont, 'Color', col.survrx);
text(x_i(1)+0.4, y_i(1)+0.4, z_i(1)+0.4, ...
     't = 0 s', 'FontSize',10, 'FontName',labelFont, 'Color',[0.25 0.25 0.30]);
text(x_i(end)-4.2, y_i(end)+0.4, z_i(end)+0.4, ...
     't = 60 s', 'FontSize',10, 'FontName',labelFont, 'Color',[0.25 0.25 0.30]);

% Baseline midpoint label
bl_mid = (Tx_Pos + RefRx_Pos)/2;
text(bl_mid(1)-2.0, bl_mid(2)+0.8, bl_mid(3)+0.5, ...
     sprintf('Baseline \\approx %.1f km', Baseline), ...
     'FontSize',10, 'FontAngle','italic', 'FontName',labelFont, ...
     'Color', col.baseline);

% -------------------- Axes cosmetics -----------------------------
xlabel('X [km]', 'FontSize',12, 'FontName',labelFont);
ylabel('Y [km]', 'FontSize',12, 'FontName',labelFont);
zlabel('Altitude Z [km]', 'FontSize',12, 'FontName',labelFont);
title('FERS Simulation Scene', 'FontSize', 14, ...
      'FontName', labelFont, 'FontWeight','bold');

xlim([-8  18]);
ylim([-2  24]);
zlim([ 0   4.5]);

view(-40, 22);

% Modern axis look — soft grid, thin frame, off-black text
ax = gca;
set(ax, 'FontSize',10.5, 'FontName',labelFont, ...
        'LineWidth', 0.9, ...
        'XColor',[0.25 0.25 0.30], 'YColor',[0.25 0.25 0.30], 'ZColor',[0.25 0.25 0.30], ...
        'GridColor',[0.55 0.58 0.65], 'GridAlpha',0.18, ...
        'MinorGridLineStyle','none', ...
        'Box','off');
% Soften the coordinate panes
try
    ax.XAxis.Color = [0.30 0.30 0.35];
    ax.YAxis.Color = [0.30 0.30 0.35];
    ax.ZAxis.Color = [0.30 0.30 0.35];
    set(ax, 'BoxStyle','back');
catch, end

% -------------------- Export -------------------------------------
outDir = fullfile(fileparts(mfilename('fullpath')), 'Figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end

% Vector SVG (best for PowerPoint)
svgPath = fullfile(outDir, 'fers_scene_scenario1.svg');
try
    exportgraphics(fig, svgPath, 'ContentType','vector', ...
                   'BackgroundColor','none');
catch  % fall back for older MATLAB
    print(fig, svgPath, '-dsvg', '-painters');
end

% High-DPI PNG fallback
pngPath = fullfile(outDir, 'fers_scene_scenario1.png');
try
    exportgraphics(fig, pngPath, 'Resolution', 300, 'BackgroundColor','none');
catch
    print(fig, pngPath, '-dpng', '-r300');
end

fprintf('wrote %s\n', svgPath);
fprintf('wrote %s\n', pngPath);
