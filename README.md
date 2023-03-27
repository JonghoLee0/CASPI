# CASPI: Collaborative Photon Processing for Active Single-Photon Imaging
This repository contains the minimum working code for our paper:

CASPI: Collaborative Photon Processing for Active Single-Photon Imaging \
Jongho Lee, Atul Ingle, JV Chacko, KW Eliceiri, Mohit Gupta \
Nature Communications 2023

<img src="https://github.com/JonghoLee0/CASPI/blob/main/teaser.png>


## Requirement
The code has been tested with MATLAB R2020b on a PC with 64-bit Windows 10.

## Demos
Run `main.m`

## Input
### Photon transient cubes
One dataset for LiDAR applicaiton and another dataset for FLIM application are provided in the folders `Data\LiDAR` and `Data\FLIM`, respectively.

### Laser source prior
FWHM values are provided in the code assuming the Gaussian laser pulses.

## Output
Recovered photon fluxes are saved in the folder `Data\LiDAR` or `Data\FLIM` according to the application. Depth maps or lifetime images are also saved in the same folder.
