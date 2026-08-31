
# Covariance Scaling Kalman Filter

A Matlab implementation of the Covariance Scaling Kalman Filter for tracking in 
the Bistatic Range and  Doppler Domain.


1.Adaptive Estimation of Measurement Noise Covariance (R): updating the measurement noise covariance matrix R adaptively based on the residual errors (ek). If the current residual error exceeds a threshold and is larger than the mean of recent errors, update the covariance matrix R.

2.Kalman Gain Calculation: You calculate the Kalman gain using the current state estimate, covariance matrices, and measurement matrix H.

4.Update State Estimate and Covariance: You update the state estimate (X) using the robustly estimated values and update the error covariance matrix (P) using the Kalman gain.

5.Reset Measurement Noise Covariance Matrix (R): Finally, reset the measurement noise covariance matrix R to predefined values
