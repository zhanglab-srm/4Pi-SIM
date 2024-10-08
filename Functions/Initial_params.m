function [dataparams,OTFparams,Reconparams] = Initial_params(pathname, filenameA, isOPLD)
    % This function will set up the initial parameters
    dataparams.pixelsizeXY = 63.73;      % rawdata xy sampling
    dataparams.pixelsizeZ = 40;          % rawdata z sampling
    dataparams.Exwave = 488;             % excitation wavelength nm
    dataparams.Emwave = 525;             % emission wavelength nm
    dataparams.NA = 1.27;                % numerical aperture 
    dataparams.RI = 1.333;               % refractive index
    dataparams.nangle = 5;               % 5 phases
    dataparams.nort = 3;                 % 3 orientations
    dataparams.spacing = [2,2,2,2,2];    % spacing for separation matrix
    dataparams.regul = 2*pi;
    dataparams.Mag = 102;                   % system magnification
    dataparams.pathname = pathname;
    dataparams.filenameA = filenameA;
    dataparams.filenameB = [filenameA(1:end-5) 'B.tif'];

    OTFparams.OTFpixelsizeZ = 40;       % OTF z sampling
    OTFparams.OTF_6B_Path='OTF\6BOTF.tif'; % Six-beam OTF
    
    Reconparams.isOPLD = isOPLD;
    Reconparams.wiener_coeff = 0.2;       % wiener constant
    Reconparams.Zpadding = 4;                % z padding
    Reconparams.bgflag = 1;                  % background noise
    Reconparams.bgname = 'Functions\background.tif';
    Reconparams.notch_switch = 1;            % notch filter
    Reconparams.notch_para_1 = 0.01;         % notch filter parameters
    Reconparams.notch_para_2 = 1.50;
    Reconparams.p_precision = 0.01;      % Subpixel precision when solving for plong
    Reconparams.p_theory = 154;          % plong theory value
    Reconparams.p_limit = 10;            % plong limit
    Reconparams.ap_limit = 150;          % apodization parameters
    Reconparams.ap_delta = 10;
    
    % Save initial parameters
    fileID = fopen([pathname 'SIM Result\Processing_Data.txt'],'w');
    fprintf(fileID,'---------------Setup_Parameters---------------\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'----------------Data_Parameters----------------\n');
    fprintf(fileID,'%1s %.2f\n','Rawdata pixelsizeXY (nm) = ',dataparams.pixelsizeXY);
    fprintf(fileID,'%1s %.2f\n','Rawdata pixelsizeZ (nm) = ',dataparams.pixelsizeZ);
    fprintf(fileID,'%1s %.2f\n','Recondata pixelsizeXY (nm) = ',dataparams.pixelsizeXY/2);
    fprintf(fileID,'%1s %.2f\n','Recondata pixelsizeZ (nm) = ',dataparams.pixelsizeZ);
    fprintf(fileID,'%1s %2d\n','Emssion Wavelength (nm) = ',dataparams.Emwave);
    fprintf(fileID,'%1s %2d\n','Excitation Wavelength (nm) = ',dataparams.Exwave);
    fprintf(fileID,'%1s %.2f\n','NA = ',dataparams.NA);
    fprintf(fileID,'%1s %.2f\n','Refractive Index = ',dataparams.RI);
    fprintf(fileID,'%1s %2d\n','Magnification = ',dataparams.Mag);    
    fprintf(fileID,'%1s %2d\n','nort = ',dataparams.nort);
    fprintf(fileID,'%1s %2d\n','nangle = ',dataparams.nangle);
    fprintf(fileID,'\n');
    fprintf(fileID,'----------------OTF_Parameters----------------\n');
    fprintf(fileID,'%1s %2d\n','isOPLD = ',isOPLD);
    fprintf(fileID,'%1s %.2f\n','OTF pixelsizeZ (nm) = ',OTFparams.OTFpixelsizeZ);
    fprintf(fileID,'\n');
    fprintf(fileID,'---------------Recon_Parameters---------------\n');
    fprintf(fileID,'%1s %.2f\n','Wiener constant = ',Reconparams.wiener_coeff);
    fprintf(fileID,'%1s %2d\n','bgflag = ',Reconparams.bgflag);
    fprintf(fileID,'%1s %2d\n','notch_switch = ',Reconparams.notch_switch);
    fprintf(fileID,'%1s %.2f\n','notch_para_1 = ',Reconparams.notch_para_1);
    fprintf(fileID,'%1s %.2f\n','notch_para_2 = ',Reconparams.notch_para_2);
    fprintf(fileID,'%1s %.2f\n','precision = ',Reconparams.p_precision);
    fprintf(fileID,'%1s %2d\n','apodization limit = ',Reconparams.ap_limit);
    fprintf(fileID,'%1s %2d\n','plong in theory = ',Reconparams.p_theory);
    fprintf(fileID,'%1s %2d\n','plong limit = ',Reconparams.p_limit);
    fprintf(fileID,'%1s %2d\n','Zpadding = ',Reconparams.Zpadding);
    fprintf(fileID,'\n');
    fprintf(fileID,'----------------Rawdata & OTF----------------\n');
    fprintf(fileID,'%1s %1s %1s\n','Rawdata_Path = ',pathname,filenameA);
    fprintf(fileID,'\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'-----------------Processing------------------\n');
    fprintf(fileID,'\n');
    fclose(fileID);
    fclose('all');
end