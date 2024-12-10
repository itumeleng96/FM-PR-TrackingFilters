# Tracking Filters for a Passive Radar System

## Author
**[Itumeleng Malemela]**  
Master's Dissertation Repository  

---

## Overview

This repository contains the code implementations of tracking filters designed for an **FM Passive Radar System**. The content is organized into:
- **Simulation code**: Includes the processing chain and all necessary (Flexible Extensible Radar Simulator) files for simulating radar scenarios.
- **Tracking Filter Algorithms**: Different Tracking filters used in the simulations.

---

## Implemented Tracking Filters

The following tracking filters and associated algorithms are implemented:

1. **Kalman Filter** 
2. **Particle Filter** 
3. **Unscented Kalman Filter**  
4. **Recursive Gauss-Newton Filter**  
5. **Covariance Scaling Techniques**

---

## Other Algorithms and Features

Additionally, this repository includes implementations of:

- **Cell-Averaging Constant False Alarm Rate (CFAR) Filter**  
- **Multi-Target Tracker**  
- **Mean Shift Clustering Algorithm**  
- **Optimization of Extensive Cancellation Algorithm**  
- **Log-Likelihood as a Statistical Benchmark** for tracking filter performance  
- **Ground Truth Calculations for FERS** using cubic waypoints  

---

## How to Run the Code

### Running a Tracking Filter
To execute a tracking filter in MATLAB, use the command:  
```matlab
runTrackingFilterPlot
```

### Running Tracking filter Evaluations
To execute the Performance evaluation script in MATLAB, use the command:  
```matlab
runEvaluateTrackingFilter
```
