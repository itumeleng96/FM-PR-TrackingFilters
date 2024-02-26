function crit = m_estimate_criterion(x, Y, X, robust_score)
    % Compute the criterion for M-estimation
    crit = 0;
    for i = 1:size(Y, 1)
        residual = Y(i, :) - X(i, :) * x;
        crit = crit + robust_score.evaluate(residual);
    end
end