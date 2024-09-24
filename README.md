# 4Pi-SIM
Matlab code for 4Pi-SIM reconstruction

# Requirements
  - Microsoft Windows 7 or newer, 64-bit
  - Matlab R2023a or newer

# How to run

Step 1: Step 1: Download the example dataset and place the 'Data' folder, which includes raw data A, data B, and OTF, in the same directory as the main code 'Main_4Pi_SIM_reconstruction.m'.

The example dataset is available at: https://doi.org/10.6084/m9.figshare.25714068.

Step 2: Open the 'Main_4Pi_SIM_reconstruction.m' code and set the initial parameter 'isOPD' to 1. Run this code. Next, select the raw data file 'ER_T86_A.tif'. The code will automatically choose the optimal OPD OTF for Wiener reconstruction. If traditional reconstruction is needed, set the initial parameter 'isOPD' to 0. In this case, the code will not estimate OPD and will reconstruct with a zero OPD OTF.

# Testing environments
  - Microsoft Windows 11 64-bit
  - Matlab R2023a
  - CPU: 12th Gen Intel(R) Core(TM) i9-12900K 3.20 GHz
  - GPU: NVIDIA GeForce RTX 3080 Ti
  - CUDA 11. 4. 136 driver

# Tested run time
225 s

# Expected output
The 'SIM result' folder will be generated, which includes an 'ER_T86_A_parameter.mat' file where the estimated parameters (wave vector, initial phase, contrast, and OPLD phase) are saved; a 'Processing_Data.txt' file where the initial and processing parameters of the code are saved; 'SpectrumWienerER_T86_A.tif', which is the spectrum of the reconstructed image, and 'WienerER_T86_A.tif', which is the final reconstructed image.

# Contact
For any questions / comments about this software, please contact zhanglab@westlake.edu.cn.

# Copyright and Software License
Copyright (c) 2024 @ Zhang Lab, Westlake University, Hangzhou, China

The package is licenced under the [GNU GPL](https://www.gnu.org/licenses/). 

