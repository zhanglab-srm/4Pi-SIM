%---------------------------------------------------------%
%%%%%%%%%%%%%%%%% 4Pi-SIM reconstruction %%%%%%%%%%%%%%%%%%
%---------------------------------------------------------%

clc
clear
close all
dbstop if error
warning off
subFolders = genpath('./Functions'); 
addpath(subFolders)
addpath('./OTF')

%% OPLD estimation
isOPLD = 1;  % 0: do not estimate OPLD; 1: estimate OPLD

%% Select raw data
currPath = fileparts(mfilename('fullpath')); 
[filenameA,pathname] = uigetfile({fullfile(currPath, '*.tif')},'choose raw data');
if isequal(filenameA,0)
    disp('User selected Cancel')
    error('User selected Cancel');
else
    disp(['User selected: ', fullfile(pathname, filenameA)])
end
mkdir([pathname 'SIM Result'])

%% Set up initial parameters
[dataparams,OTFparams,Reconparams] = Initial_params(pathname, filenameA, isOPLD);
diary([pathname 'SIM Result\Processing_Data.txt'])
clear pathname filenameA isOPLD subFolders currPath 

%% Perform SIM reconstruction
[dataparams, OTFparams, Reconparams] = EstimateParameters(dataparams, OTFparams, Reconparams);

Wiener_deconv_4PS(dataparams, OTFparams, Reconparams);

diary off
