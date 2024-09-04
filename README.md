# Tracking Filters For a Passive Radar System

### This Repository will consist Code implementations for Tracking Filters for a FM Passive Radar System.<br>

#### The Repository is divided into code processing chain for Simulations and all the required files and Real world data folder for the Radar Processing Chain For real-data

#### The Following Tracking filters and code have been implemented:
<li>1.Kalman Filter and Adaptive Kalman Filter</li>
<li>2.Particle Filter  and Adaptive Particle Filter</li>
<li>3.Unscented Kalman and Adaptive Unscented Kalman Filter </l>
<li>4.Recursive Gauss Newton Filter and Adaptive Recursive Gauss Newton Filter</li>


#### Other algorithms implemented
<li>Constant False Alarm Rate Filter </li>
<li>Multi-target Tracker</li>
<li>MeanshiftCluster</li>
<li>Log-likelihood as a statistical Benchmark for the Performance of Tracking Filters </li>
<li>Ground truth calculations for FERS for cubic waypoints</li>


### Planned Activity
<li>Tracker optimization based on Data from FM based Peralex PR deployed at the SKA.</li>
<li>Target based tracker optimization in such a way that the tracking filter parameters can be tuned based on the target of interest</li>
<li>An investigation on the possibility of running multiple concurrent trackers if there are multiple targets of interest that have different target characteristics.</li>
<li>A comprehensive report detailing the optimal parameters, advantages and disadvantages of each filter and the various filter parameters</li>
<li>Accept real-time streamed data into the Multi-Target Tracking Filter if time permits.</li>

## How to run the tracking Filter in the Matlab Terminal
<li> runTrackingFilter </li>
<li> Enter the filterType:  1,2,3,4,5,6,7 for Kalman,Particle,Adaptive Kalman,UKF,AUKF,RGNF,ARGNF filter </li>


## Filter Evaluations with Log-Likelihood
<li> runEvaluateFilters </li>
