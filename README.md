# gpr_ms_thesis
MS Thesis scripts from my time at University of Rhode Island

## Background
Gaussian Process Regression (GPR) was developed in 1995 as a machine learning method to solve the regression problem. GPR is extremely accurate and it provides uncertainty estimates of each prediction point which makes it more valuable than comparable methods. GPR suffers from O(N^3) time complexity, and thus researchers have posed numerous ways to approximate the core equations to speed things up. Kris Krasnosky did his PhD at U Rhode Island (URI) on the development of an online GPR solution that utilized a GPU and other computations speed-ups. Being in the Oceanography field, Krasnosky focused on making 3D maps of the seafloor.

## Problem
The problem is twofold.
1. GPR is still a heavy computational lift even with a GPU under Krasnosky's approach. Plus, we're dealing with a massive dataset (multi-beam sonar data is dense). Can we make GPR faster?
2. Not all sonar data is created equal. The physics that govern underwater acoustic waves cause the outer beams to be less accurate than inner beams. Can we include a physics model to more accurately represent sonar data in GPR?

## Approach
For each problem there is an approach to solving it.
1. There are numerous approximate methods that exist to reduce the input data to GPR, which reduces the compute time. I will test a variety of these downsampling methods and see which reduces the compute time the most while holding accuracy high (when compared to standard GPR).
2. There are papers that propose methods for quantifying the uncertainty associated with each sonar beam in a multibeam array. This should serve as an input to the GPR algorithm, potentially replacing the sensor variance (a GPR hyperparameter) with something that captures the range of uncertainty.

## Generalizing
This Ocean Engineering MS Thesis is built and tuned for multi-beam sonar data. The scripts in this repository are for preliminary testing and experimentation. The actual full scale implementation (C++, ROS) is kept within the Roman Lab at URI.

With that being said, the approaches used here can be generalized to other research in the environmental variable modeling domain. 

## Results
Results will be presented in my MS thesis and in a publication that I will link at a later date (~June 2023).

## Acknowledgements
Shoutout to the original researchers who developed GPR, to Kris Krasnosky for bringing it into the modern era and to URI, and to Chris Roman for advising me on this project. 
