# Automatic Gait Events Detection with IMU: An in-Depth Look at the MP-DTW Algorithm

## General informations 

IPOL link: 
Date: 
Author mail: cyril.voisard@gmail.com

Please cite these articles whenever you want to make a reference to this demo:
- [Voisard et al, Automatic gait events detection with inertial measurement units: healthy subjects and moderate to severe impaired patients, Journal of NeuroEngineering and Rehabilitation, 2024, https://doi.org/10.1186/s12984-024-01405-x]
- [Voisard et al., Automatic Gait Events Detection with IMU: An in-Depth Look at the MP-DTW Algorithm, IPOL, 2025, XXXXX]

## Context

This repository provides the code for the MP-DTW algorithm, a tool for automated gait event detection designed for clinical tests with IMUs recordings. 

The first stage of scientific validation of this tool was published in 2024 in the Journal of NeuroEngineering and Rehabilitation [Voisard et al, Automatic gait events detection with inertial measurement units: healthy subjects and moderate to severe impaired patients, Journal of NeuroEngineering and Rehabilitation, 2024, https://doi.org/10.1186/s12984-024-01405-x]. 
This article explains all the details about the algorithm to allow reproducibility.

An article detailing each of the mathematical formulas used and the associated lines of code has also been published in IPOL : 
[Voisard et al., Automatic Gait Events Detection with IMU: An in-Depth Look at the MP-DTW Algorithm, IPOL, 2025, XXXXX].


## Demo and accessibility 

The demo applying the segmentation to data supplied in the correct format is available here: https://www.ipol.im/pub/art/XXXX.
The format of the input data and the output figures and tables are provided in the attached article. 


## Examples of use

3 datasets are provided in the folder "example_datasets" consisting of 2 gait signals collected using two inertial measurement units (MTw Awinda XSens, on each dorsal part of the feet). Two trials (noted 2 x 10m) was obtained from a sample of subjects who followed the described protocol: standing still, walking 10 meters, turning around, walking back, and stopping. The third trials was a test with 46 3-meters straight line, with a u-turn between each. 
The data output by the algorithm is provided in the same files. 
