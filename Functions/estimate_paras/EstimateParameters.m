function [dataparams, OTFparams, Reconparams] = EstimateParameters(dataparams,OTFparams,Reconparams)
    % This function will estimate plong, initial phases, contrast and OPLD.

    % load SIM Raw data (XYAPZ)
    disp('Loading SIM Rawdata...')
    endframe = numel(imfinfo([dataparams.pathname dataparams.filenameA]))/ ...
        (dataparams.nangle*dataparams.nort);
    SIM_rawA = imreadstack_TIRF([dataparams.pathname dataparams.filenameA], 1,...
        endframe*dataparams.nangle*dataparams.nort);
    SIM_rawB = imreadstack_TIRF([dataparams.pathname dataparams.filenameB], 1,...
        endframe*dataparams.nangle*dataparams.nort);
    [sizex,sizey,sizez] = size(SIM_rawA);
    nangle = dataparams.nangle;
    nort = dataparams.nort;
    numz = sizez/(nangle*nort);
    Zmedian = floor(numz/2)+1;
    max_xy = max(sizex,sizey);

    dataparams.sizex = sizex;
    dataparams.sizey = sizey;
    dataparams.numz = numz;       

    disp('Estimating SIM parameters...')
    Pindex = 1;
    Phalf_index = 2;
    if Reconparams.isOPLD == 1 % estimate OPLD
        % registration for SIM rawdataA and rawdataB
        [SIM_raw6B,~] = Run_4ps_registration(SIM_rawA, SIM_rawB);
        clear SIM_rawA SIM_rawB
        
        % resampling and normalization
        H = imreadstack(OTFparams.OTF_6B_Path);
        H_Resample = Resample_OTF(H, OTFparams.OTFpixelsizeZ, dataparams.pixelsizeZ, numz, max_xy);
           
        FitBand.H0_para = H_Resample(:,:,Zmedian);                            
        FitBand.H2_para = H_Resample(:,:,numz*Pindex+Zmedian);
        FitBand.H1_para = H_Resample(:,:,numz*Phalf_index+Zmedian);   
        pz = dataparams.pixelsizeZ * numz * dataparams.RI / dataparams.Exwave;
        pz = round(pz);
        inteval = round(pz/2);
        pz_opld = round(Zmedian-3/2*inteval);    % the layer includes OPLD
        FitBand.H0_opld = H_Resample(:,:,pz_opld);                             
        FitBand.H1_opld = H_Resample(:,:,numz*Phalf_index + pz_opld);
        clear H H_Resample pz inteval
    
        % Parameters estimation
        Reconparams.coodinatex = zeros(9,2);
        Reconparams.coodinatey = zeros(9,2);
        Reconparams.initial_phase = zeros(nort,2);
        Reconparams.contrast = zeros(nort,2);
        phase = zeros(nort,1);
        Reconparams.opld = zeros(nort,1);
        for orti=1:nort
            % separate bands in each direction
            [sep_im] = Separate_SIMbands(SIM_raw6B(:,:,1:endframe*nangle*nort), ...
                nangle, nort, dataparams.spacing, dataparams.regul, numz, orti); 
            
            % fitting plong, initial phases and contrast
            FitBand.band0_para = sep_im(:,:,Zmedian,3);
            FitBand.band2_para = sep_im(:,:,Zmedian,Pindex);
            FitBand.band1_para = sep_im(:,:,Zmedian,Phalf_index);
            sep_im = sep_im./(abs(sep_im)+eps);
            FitBand.band0_p = sep_im(:,:,Zmedian,3);
            FitBand.band2_p = sep_im(:,:,Zmedian,Pindex);
            FitBand.band0_a = sep_im(:,:,pz_opld,3);
            FitBand.band1_a = sep_im(:,:,pz_opld,Phalf_index);
            clear sep_im 
            
            [tmpx,tmpy,tmpphase,tmpc] = Fitting_SIMpara(FitBand, dataparams, Reconparams, orti);

            Reconparams.coodinatex((orti-1)*3+1:(orti-1)*3+3,:) = tmpx((orti-1)*3+1:(orti-1)*3+3,:);
            Reconparams.coodinatey((orti-1)*3+1:(orti-1)*3+3,:) = tmpy((orti-1)*3+1:(orti-1)*3+3,:);
            Reconparams.initial_phase(orti,:) = tmpphase;
            Reconparams.contrast(orti,:) = tmpc;
            plong = ((Reconparams.coodinatex((orti-1)*3+1:(orti-1)*3+3,:)-Reconparams.coodinatex(1,1)).^2+ ...
                (Reconparams.coodinatey((orti-1)*3+1:(orti-1)*3+3,:)-Reconparams.coodinatey(1,1)).^2).^0.5;

            disp('------------------------------------')
            disp([' From ori' num2str(orti) ':       order2   order1']);
            disp([' plong:           ' num2str(plong(2,1),'%.2f') '   ' num2str(plong(2,2),'%.2f')]);
            disp([' contrast:        ' num2str(Reconparams.contrast(orti,1),'%.2f') '     ' ...
                num2str(Reconparams.contrast(orti,2),'%.2f')]);
            disp([' initial phase:   ' num2str(Reconparams.initial_phase(orti,1),'%.2f') '     ' ...
                num2str(Reconparams.initial_phase(orti,2),'%.2f')]);
    
            % fitting opld  
            [tmpphase] = Fitting_OPLD(FitBand, dataparams, Reconparams, orti);

            phase(orti,:) = tmpphase(1,2);
            Reconparams.opld(orti,:) = phase(orti,:)*dataparams.Exwave/dataparams.RI/2/pi;
            savedata = rmfield(Reconparams, setdiff(fieldnames(Reconparams), {'coodinatex','coodinatey', ...
                'initial_phase','contrast','opld'}));
            save([dataparams.pathname 'SIM Result/' dataparams.filenameA(1:end-4) ...
                '_parameter.mat'], '-struct', 'savedata');
        end
        disp('------------------------------------')
        OPLD_total = sum(Reconparams.opld,"all")./nort;
        Reconparams.OPLD_total = OPLD_total;
        disp([' OPLD : ' num2str(OPLD_total) ' nm']);
        disp('------------------------------------')
        
        clear SIM_raw6B endframe nangle nort pz_opld plong phase sizex sizey sizez tmpc ...
         tmpphase tmpx tmpy Zmedian Pindex Phalf_index orti

        % choose OTF with correct OPLD: We provide two representative OTFs selected from the 4Pi-SIM OTF library.
        if 20 <= OPLD_total && OPLD_total < 30  % OPLD = 20 nm
            H_real = imreadstack('OTF\4PSOTF_OPLD_20_real.tif');
            H_imag = imreadstack('OTF\4PSOTF_OPLD_20_imag.tif');
        elseif 180 <= OPLD_total && OPLD_total < 190  % OPLD = 180 nm
            H_real = imreadstack('OTF\4PSOTF_OPLD_180_real.tif');
            H_imag = imreadstack('OTF\4PSOTF_OPLD_180_imag.tif');
        end
        H_complex = complex(H_real,H_imag);
        clear H_real H_imag
        OTFparams.H_Recon = Resample_OTF(H_complex, OTFparams.OTFpixelsizeZ, dataparams.pixelsizeZ, numz, max_xy);
        clear H_complex
        
    else  % do not estimate OPLD
        H_real = imreadstack('OTF\4PSOTF_OPLD_0_real.tif');
        H_imag = imreadstack('OTF\4PSOTF_OPLD_0_imag.tif');
        H_complex = complex(H_real,H_imag);
        % resampling and normalization
        OTFparams.H_Recon = Resample_OTF(H_complex, OTFparams.OTFpixelsizeZ, dataparams.pixelsizeZ, numz, max_xy);
            
        FitBand.H0_para = OTFparams.H_Recon(:,:,Zmedian);                             
        FitBand.H2_para = OTFparams.H_Recon(:,:,numz*(Pindex)+Zmedian);
        FitBand.H1_para = OTFparams.H_Recon(:,:,numz*Phalf_index+Zmedian); 
        clear H_complex H_real H_imag 

        % parameters estimation
        Reconparams.coodinatex = zeros(9,2);
        Reconparams.coodinatey = zeros(9,2);
        Reconparams.initial_phase = zeros(nort,2);
        Reconparams.contrast = zeros(nort,2);
        for orti=1:nort
            % separate bands in each direction
            [sep_im] = Separate_SIMbands(SIM_rawA(:,:,1:endframe*nangle*nort), ...
                nangle, nort, dataparams.spacing, dataparams.regul, numz, orti); 
  
            % fitting plong, initial phases and contrast
            FitBand.band0_para = sep_im(:,:,Zmedian,3);            
            FitBand.band2_para = sep_im(:,:,Zmedian,Pindex);          
            FitBand.band1_para = sep_im(:,:,Zmedian,Phalf_index); 
            sep_im = sep_im./(abs(sep_im)+eps);
            FitBand.band0_p = sep_im(:,:,Zmedian,3);
            FitBand.band2_p = sep_im(:,:,Zmedian,Pindex);
            clear sep_im
    
            [tmpx,tmpy,tmpphase,tmpc] = Fitting_SIMpara(FitBand, dataparams, Reconparams, orti);
            
            Reconparams.coodinatex((orti-1)*3+1:(orti-1)*3+3,:) = tmpx((orti-1)*3+1:(orti-1)*3+3,:);
            Reconparams.coodinatey((orti-1)*3+1:(orti-1)*3+3,:) = tmpy((orti-1)*3+1:(orti-1)*3+3,:);
            Reconparams.initial_phase(orti,:) = tmpphase;
            Reconparams.contrast(orti,:) = tmpc;
            plong = ((Reconparams.coodinatex((orti-1)*3+1:(orti-1)*3+3,:)-Reconparams.coodinatex(1,1)).^2+ ...
                (Reconparams.coodinatey((orti-1)*3+1:(orti-1)*3+3,:)-Reconparams.coodinatey(1,1)).^2).^0.5;

            disp('------------------------------------')
            disp([' From ori' num2str(orti) ':       order2   order1']);
            disp([' plong:           ' num2str(plong(2,1),'%.2f') '   ' num2str(plong(2,2),'%.2f')]);
            disp([' contrast:        ' num2str(Reconparams.contrast(orti,1),'%.2f') '     ' ...
                num2str(Reconparams.contrast(orti,2),'%.2f')]);
            disp([' initial phase:   ' num2str(Reconparams.initial_phase(orti,1),'%.2f') '     ' ...
                num2str(Reconparams.initial_phase(orti,2),'%.2f')]);
        end
        disp('------------------------------------')
        clear SIM_raw6B endframe nangle nort plong phase sizex sizey sizez tmpc ...
            tmpphase tmpx tmpy Zmedian Pindex Phalf_index orti

        savedata = rmfield(Reconparams, setdiff(fieldnames(Reconparams), {'coodinatex','coodinatey', ...
            'initial_phase','contrast'}));
        save([dataparams.pathname 'SIM Result/' dataparams.filenameA(1:end-4) ...
            '_parameter.mat'], '-struct', 'savedata');
    end
    
    % calculate apodization parameter
    Reconparams.ap_range = dataparams.pixelsizeXY*100/dataparams.Mag*max_xy*2*dataparams.NA/ ...
        dataparams.Emwave+Reconparams.ap_delta;
end