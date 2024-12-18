# Tracking Filters For a Passive Radar System

### This Repository will consist Code implementations for Tracking Filters for a FM Passive Radar System.<br>

#### The Repository is divided into code processing chain for Simulations and all the required files and Real world data folder for the Radar Processing Chain For real-data

#### The Following Tracking filters and code have been implemented:
<li>1.Kalman Filter and Adaptive Kalman Filter</li>
<li>2.Particle Filter  and Adaptive Particle Filter</li>
<li>3.Unscented Kalman and Adaptive Unscented Kalman Filter </li>
<li>4.Recursive Gauss Newton Filter and Adaptive Recursive Gauss Newton Filter</li>
<li>5.Covariance Scaling Techniques</li>

#### Other algorithms implemented
- Constant False Alarm Rate Filter  

  ![CFAR Results](/Simulation_results/CFAR/cacfar_pfa-8.svg)

- Multi-target Tracker

    <video controls width="485">
        <source src="/Simulation_results/mttOutput.mp4" type="video/mp4">
        Your browser does not support the video tag.
    </video>

- MeanshiftCluster 
  ![CFAR Centroid Results](/Simulation_results/CFAR/centroidCFAR.png)

- Log-likelihood as a statistical Benchmark for the Performance of Tracking Filters

   ![Log Likelihood](/Simulation_results/FERS_scenarios/360_range_ll_1.svg)

- Ground truth calculations for FERS for cubic waypoints

   ![GT](/Simulation_results/FERS_scenarios/3D_360.svg)
   ![GT](/Simulation_results/FERS_scenarios/rangeDoppler360.svg)

## How to run the tracking Filter in the Matlab Terminal
<li> runTrackingFilterPlot </li>

## Filter Evaluations with Log-Likelihood
<li> runEvaluateFilters </li>
