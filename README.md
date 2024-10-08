# 4Pi-SIM
Matlab code for 4Pi-SIM reconstruction

# Requirements
  - Microsoft Windows 7 or newer, 64-bit
  - Matlab R2023a or newer

# How to run

Step 1: Download the example dataset and place the raw data and the 'OTF' folder in the same directory as the main code, 'Main_4Pi_SIM_reconstruction.m'.

The example dataset is available at: https://doi.org/10.6084/m9.figshare.25714068.

Step 2: Open 'Main_4Pi_SIM_reconstruction.m', set the initial parameter 'isOPLD=1', and run the code. The code will automatically select the OTF with the optimal OPLD for Wiener reconstruction when the raw data is chosen. For standard reconstruction, set the initial parameter to 'isOPLD=0' to bypass OPLD estimation and use an OTF with zero OPLD for data reconstruction.

# Testing environments
  - Microsoft Windows 11 64-bit
  - Matlab R2023a
  - CPU: 12th Gen Intel(R) Core(TM) i9-12900K 3.20 GHz
  - GPU: NVIDIA GeForce RTX 3080 Ti
  - CUDA 11. 4. 136 driver

# Tested run time
197.6 s

# Expected output
  - Fig5a_ER_4Pi_SIM_A_parameter.mat: Estimated parameters (wave vector, initial phase, contrast, and OPLD phase)
  - Processing_Data.txt: Initial settings and processing parameters
  - SpectrumWienerFig5a_ER_4Pi_SIM_A.tif: Spectrum of the reconstructed image
  - WienerFig5a_ER_4Pi_SIM_A.tif: Final reconstructed image 


# Contact
For any questions/comments about this code, please contact zhanglab@westlake.edu.cn.

# Copyright and Software License
Copyright (c) 2024 @ Zhang Lab, Westlake University, Hangzhou, China

The package is licensed under the [GNU GPL](https://www.gnu.org/licenses/). 



