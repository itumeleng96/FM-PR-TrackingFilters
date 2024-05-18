function ll = logLikelihood(mean, S ,sample)
    % Inputs:
    %   - sample: Array of observed  values from the measured track.
    %   - mean: Mean of the distribution under consideration from the predicted track.
    %   - S: Innovation Covariance value for the parameter (scalar).
    % 
    % Output:
    %   - ll: Log-likelihood for the parameter.
    
    disp(mean);
    disp(S);
    disp(sample);

   
    % Calculate the log-likelihood for the  parameter
    
    % Log-likelihood for a NORMAL distribution is given by:
    ll = -0.5 * log(2 * pi * S^2) - ((sample - mean)^2) / (2 * S^2);
    %ll = log(normpdf(sample,mean,S));
end