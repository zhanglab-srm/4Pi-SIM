%----------------------------------------------------------
% 4Pi-SIM reconstruction
%----------------------------------------------------------
clc
clear
close all
dbstop if error
warning off
addpath('./Functions')
%% Initial parameters
isOPLD = 1;  % 0: dont estimate OPLD; 1: estimate OPLD
pixelsize = 63.73;        % rawdata xy sampling
pixelsizeZ = 40;          % rawdata z sampling
Beads_pixesizeZ = 40;     % OTF z sampling
weilac = 0.2;             % wiener coefficient
wavelengh = 488;          % excitation wavelength nm
EmissionWavelength = 525; % emission wavelength nm
NA = 1.27;
Pzrefmed = 1.333;
nangle = 5;          % 5 phases
nort = 3;            % 3 orientations
spjg = [2,2,2,2,2];  % 1:1:1:1:1
regul = 2*pi;
Zpadding = 4;
Mode = '3D-SIM-3';  % imaging mode
M = 102;            % System magnification

% Six-beam OTF
OTF_6B_Path='Data\OTF\Assemble-6BOTF-abs.tif';

bgflag = 1;            % background noise
bgname = 'Functions\background.tif';
notch_swith = 1;       % notch filter
notch_para = 10;       % notch filter parameter
notch_para_1 = 0.01;
notch_para_2 = 1.50;
jindu = 0.01;      % Subpixel accuracy when solving for plong in cross-correlation
pg = 154;          % plong theory value
fanwei = 10;       % plong limit
fc_para = 150;     % apodization parameters
fc_delta = 10;
I2M_para = 5;
isEqualization = 1; % equalise the 15 z-stacks
isNormalizOTF = 1;  % normalise the OTF after resampling
isSaveWF = 0;       % WF images

%% Select data
currPath = fileparts(mfilename('fullpath')); % get current path
p = cd(currPath);
[filenameA,pathname]=uigetfile({'D:\210705_WCY_WY\*.tif'},'choose raw data');
if isequal(filenameA,0)
    disp('User selected Cancel')
    error('User selected Cancel');
else
    disp(['User selected: ', fullfile(pathname, filenameA)])
end

filename_all = dir([pathname '*.tif']);
number = length(filename_all);
for filei = number:-1:1
    if contains(filename_all(filei).name,'B')
        filename_all(filei) = [];
    end
end
number = length(filename_all);
mkdir([pathname 'SIM Result'])

%% Save initial parameters
diary([pathname 'SIM Result\Processing_Data.txt'])
fileID = fopen([pathname 'SIM Result\Processing_Data.txt'],'w');
fprintf(fileID,'---------------Setup_Parameters---------------\n');
fprintf(fileID,'%1s %.4f\n','pixelsizeXY = ',pixelsize);
fprintf(fileID,'%1s %.4f\n','pixelsizeZ = ',pixelsizeZ);
fprintf(fileID,'%1s %.4f\n','weilac = ',weilac);
fprintf(fileID,'%1s %2d\n','Ex_wavelength = ',wavelengh);
fprintf(fileID,'%1s %2d\n','EmWavelength = ',EmissionWavelength);
fprintf(fileID,'%1s %.2f\n','NA = ',NA);
fprintf(fileID,'%1s %.2f\n','Pzrefmed = ',Pzrefmed);
fprintf(fileID,'%1s %2d\n','bgflag = ',bgflag);
fprintf(fileID,'%1s %2d\n','notch_swith = ',notch_swith);
fprintf(fileID,'%1s %.2f\n','notch_para = ',notch_para);
fprintf(fileID,'%1s %.2f\n','notch_para_1 = ',notch_para_1);
fprintf(fileID,'%1s %.2f\n','notch_para_2 = ',notch_para_2);
fprintf(fileID,'%1s %.2f\n','jindu = ',jindu);
fprintf(fileID,'%1s %.2f\n','Beads_pixesizeZ = ',Beads_pixesizeZ);
fprintf(fileID,'%1s %.2f\n','fc = ',fc_para);
fprintf(fileID,'%1s %.2f\n','pg = ',pg);
fprintf(fileID,'%1s %.2f\n','fanwei = ',fanwei);
fprintf(fileID,'%1s %2d\n','nangle = ',nangle);
fprintf(fileID,'%1s %2d\n','nort = ',nort);
fprintf(fileID,'%1s %1s %1s %1s %1s %2d\n','spjg = ',spjg);
fprintf(fileID,'%1s %2d\n','Zpadding = ',Zpadding);
fprintf(fileID,'%1s %2d\n','Magnification = ',M);
fprintf(fileID,'%1s %2d\n','I2M_para = ',I2M_para);
fprintf(fileID,'%1s %2d\n','isEqualization = ',isEqualization);
fprintf(fileID,'%1s %2d\n','isNormalizOTF = ',isNormalizOTF);
fprintf(fileID,'%1s %2d\n','isSaveWF = ',isSaveWF);
fprintf(fileID,'\n');
fprintf(fileID,'---------------Rawdata & OTF---------------\n');
fprintf(fileID,'%1s %1s %1s\n','Rawdata_Path = ',pathname,filenameA);
fclose(fileID);
fclose('all');

%% SIM reconstruction
for filei=1:number
    filenameA = filename_all(filei).name;
    filenameB = [filenameA(1:end-5) 'B.tif'];
    SIM_reconstruction(pixelsize,pixelsizeZ,weilac,wavelengh,EmissionWavelength,NA,Pzrefmed, ...
        bgflag,notch_swith,notch_para,notch_para_1,notch_para_2,jindu,bgname,Beads_pixesizeZ,fc_para,fc_delta,pg,fanwei,nangle, ...
        nort,spjg,regul,Zpadding,M,I2M_para,isEqualization,isNormalizOTF,isSaveWF,pathname,filenameA,filenameB,OTF_6B_Path,isOPLD);
end

diary off