function ll = logLikelihood(mean, R ,sample)
    % Inputs:
    %   - sample: Array of observed  values from the predicted track.
    %   - mean: Mean of the distribution under consideration from the true track.
    %   - R: Covariance value for the parameter (scalar).
    % 
    % Output:
    %   - ll: Log-likelihood for the parameter.

   
    % Calculate the log-likelihood for the  parameter
    
    % Log-likelihood for a NORMAL distribution is given by:
    ll = log(normpdf(sample,mean,R));

end