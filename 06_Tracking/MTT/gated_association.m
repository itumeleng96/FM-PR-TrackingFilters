function [z, associated_idx, eta_all] = gated_association(clusterCentroids, x_pred, H, S, gamma)
% GATED_ASSOCIATION  Global-Nearest-Neighbour gating for mean-shift centroids.
%
%   Picks the centroid closest (in Mahalanobis sense) to the predicted
%   measurement that also passes a chi-squared gate. Returns z=[NaN;NaN]
%   if no centroid passes (caller should coast on prediction).
%
%   Inputs
%     clusterCentroids : 2xM matrix, [range_km; doppler_Hz] per column
%     x_pred           : predicted state vector (from filter.predict())
%     H, S             : measurement matrix, innovation covariance
%     gamma            : chi-squared gate threshold (2 dof, e.g. 9.21 @ 99%)
%
%   Outputs
%     z              : 2x1 measurement or [NaN;NaN] if no gated pick
%     associated_idx : column index into clusterCentroids of the pick (0 if none)
%     eta_all        : 1xM Mahalanobis distances (for diagnostics)

    if isempty(clusterCentroids)
        z = [NaN; NaN]; associated_idx = 0; eta_all = [];
        return;
    end

    z_pred = H * x_pred;
    M      = size(clusterCentroids, 2);
    eta_all = nan(1, M);
    for i = 1:M
        d = clusterCentroids(:, i) - z_pred;
        eta_all(i) = d.' * (S \ d);
    end

    [eta_min, idx] = min(eta_all);
    if eta_min <= gamma
        z = clusterCentroids(:, idx);
        associated_idx = idx;
    else
        z = [NaN; NaN];
        associated_idx = 0;
    end
end
