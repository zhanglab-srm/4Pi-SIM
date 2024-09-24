# 4Pi-SIM
Matlab code for 4Pi-SIM reconstruction

Requirements
Microsoft Windows 7 or newer, 64-bit
Matlab R2023a or newer

How to run

Step 1: Place the ‘Data’ folder (including raw dataA, dataB and OTF) in the same folder as the main code ‘Main_4Pi_SIM_reconstruction.m’.

Download links for ‘Data’ folder: https://doi.org/10.6084/m9.figshare.25714068.

Step 2: Open ‘Main_4Pi_SIM_reconstruction.m’ code. Let the initial parameter ‘isOPD=1’, and run this code then select the raw data ‘ER_T86_A.tif’, the code will automatically select the optimal OPLD OTF for Wiener reconstruction. If ordinary reconstruction is needed, you can let the initial parameter ‘isOPD=0’, then this code won’t estimate OPLD, just reconstruct with zero OPLD OTF.

Testing environments
Microsoft Windows 11 64-bit
Matlab R2023a
CPU: 12th Gen Intel(R) Core(TM) i9-12900K 3.20 GHz
GPU: NVIDIA GeForce RTX 3080 Ti
CUDA 11. 4. 136 driver

Tested run time
225 s

Expected output
One ‘SIM result’ folder includes: a ‘ER_T86_A_parameter.mat’ file where estimated parameters (wave vector, initial phase, contrast and OPLD phase) are saved; a ‘Processing_Data.txt’ file where code initial and processing parameters are saved; ‘SpectrumWienerER_T86_A.tif’ is reconstructed image’s spectrum and ‘WienerER_T86_A.tif is final reconstructed image.
