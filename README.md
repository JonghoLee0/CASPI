# CASPI: Collaborative Photon Processing for Active Single-Photon Imaging
This repository contains the minimum working code for our paper:

**CASPI: Collaborative Photon Processing for Active Single-Photon Imaging** \
Jongho Lee, Atul Ingle, JV Chacko, KW Eliceiri, Mohit Gupta \
Nature Communications 2023

\
\
<img src="https://github.com/JonghoLee0/CASPI/blob/main/teaser.png" width="800">
<br />
<br />
<br />


## Requirements
* The code has been tested with MATLAB R2020b on a PC with 64-bit Windows 10.


## Demos
* Run `main.m`


## Input
### Noisy and distorted raw photon transient cubes
* Example datasets for LiDAR applicaiton and FLIM application are provided in the folders `Data\LiDAR` and `Data\FLIM`, respectively.


### Laser source priors
* FWHM values of the Gaussian laser pulses are provided in the code for the provided datasets.
* For the datasets using non-Gaussian laser pulses, provide the number of frequency bins of the laser pulse for the parameter `N_sig_f` manually.


### Intensity images (optional)
* Intensity images (`i_map_set`) for the LiDAR dataset are provided.


## Output
* Recovered photon transient cubes are saved in the folder `Data\LiDAR` or `Data\FLIM` according to the application.
* Depth maps or lifetime images are saved in the same folder.
