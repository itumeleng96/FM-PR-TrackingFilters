function crlb = cramerRLB(N, R)
 % Inputs:
    %   - trueTrack: Array of true track values (non-log-transformed).
    %   - R: Covariance value for the parameter (scalar).
    % 
    % Output:
    %   - crlb: Cramér-Rao Lower Bound for the parameter (mean).

    % Log-transform the trueTrack values to prepare for log-likelihood calculations

    % Get the number of observations in the log_trueTrack
    num_observations = N;

    % Calculate the Fisher Information for the mean parameter of the log-normal distribution
    % The Fisher Information for the mean of a log-normal distribution with covariance R
    % is given by: I_mean = num_observations / (2 * R).

    % Invert the Fisher Information to obtain the Cramér-Rao Lower Bound (CRLB)
    
    crlb = R/num_observations;
end