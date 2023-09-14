# Tracking Filters For a Passive Radar System

### This Repository will consist Code implementations for Tracking Filters for a FM Passive Radar System.<br>


#### The Following Tracking filters and code have been implemented:
<li>1.Kalman Filter</li>
<li>2.Particle Filter  </li>
<li>3.Unscented Kalman Filter -In progress</li>
<li>4.Recursive Gauss Newton Filter - Still to be implemented</li>



#### Other algorithms implemented
<li>Constant False Alarm Rate Filter </li>
<li>Multi-target Tracker</li>
<li>MeanshiftCluster</li>
<li>Log-likelihood as a statistical Benchmark for the Performance of Tracking Filters </li>

### Planned Activity
<li>Implement other tracking filters </li>
<li>Tracker optimization based on Data from FM based Peralex PR deployed at the SKA.</li>
<li>Target based tracker optimization in such a way that the tracking filter parameters can be tuned based on the target of interest</li>
<li>An investigation on the possibility of running multiple concurrent trackers if there are multiple targets of interest that have different target characteristics.</li>
<li>A comprehensive report detailing the optimal parameters, advantages and disadvantages of each filter and the various filter parameters</li>
<li>Accept real-time streamed data into the Multi-Target Tracking Filter if time permits.</li>

## To Do's
<li>Ground Truth Data for Multiple waypoints</li>
<li>Validate Filter results</li>
<li>Implement UKF</li>
<li>Test MTT</li>

## How to run the tracking Filter in the Matlab Terminal
<li> runTrackingFilter </li>
<li> Enter the filterType:  1 for Kalman and 2 for Particle filter </li>


## Filter Evaluations
<li> runEvaluateFilters </li>

![Screenshot 1](https://github.com/itumeleng96/trackingFilters/blob/main/FilterEvaluationsLogLikelihood.png)
![Screenshot 2](https://github.com/itumeleng96/trackingFilters/blob/main/FilterEvaluationErrors.png)
