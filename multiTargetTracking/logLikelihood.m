function ll = logLikelihood(mean, R ,sample)
    % Inputs:
    %   - sample: Array of observed  values from the predicted track.
    %   - mean: Mean of the distribution under consideration from the true track.
    %   - R: Covariance value for the parameter (scalar).
    % 
    % Output:
    %   - ll: Log-likelihood for the parameter.

    % Get the number of observations in the range_sample
    num_observations = length(sample);
    
    % Calculate the log-likelihood for the  parameter
    
    % Log-likelihood for a NORMAL distribution is given by:
    % log_likelihood = -((N / 2) * log(2*pi*sigma^2)) - sum((x - mu)^2) / (2*sigma^2)
    %ll = ((num_observations/ 2) * log(2*pi*R.^2) - sum((sample - mean)^2) / (2*R.^2));


    % Log-likelihood for  a LOG-NORMAL  distribution is given by:
    % log_likelihood = -((N / 2) * log(2*pi*sigma^2)) - sum((log(x) - mu)^2) / (2*sigma^2)  
    ll = ((num_observations / 2) * log(2*pi*R.^2) - sum((log(sample) - mean).^2) / (2*R.^2));

end