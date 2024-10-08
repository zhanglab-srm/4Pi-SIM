function [] = SIM_reconstruction(pixelsize,pixelsizeZ,weilac,wavelengh,EmissionWavelength,NA,Pzrefmed, ...
    bgflag,notch_swith,~,notch_para_1,notch_para_2,jindu,bgname,Beads_pixesizeZ,fc_para,fc_delta,pg,fanwei,nangle, ...
    nort,spjg,regul,Zpadding,M,I2M_para,isEqualization,isNormalizOTF,isSaveWF,pathname,filenameA,filenameB,OTF_6B_Path,isOPLD)

%% Read data
endframe = numel(imfinfo([pathname filenameA]))/(nangle*nort);
SIM_rawA = imreadstack_TIRF([pathname filenameA],1,endframe*nangle*nort);
SIM_rawB = imreadstack_TIRF([pathname filenameB],1,endframe*nangle*nort);
[sizex,sizey,sizez] = size(SIM_rawA);
numz = sizez/(nangle*nort);
Zmedian = floor(numz/2)+1;
max_xy = max(sizex,sizey);

%% Save WF image
if isSaveWF == 1
    temp = reshape(SIM_rawA,max_xy,max_xy,nangle*nort,numz);
    WF_image = squeeze(mean(temp,3));
    imwritestack(WF_image,[pathname 'SIM Result/Widefield-' filenameA(1:end-4) '.tif'])
    WF_image_fft = fftshift(fftn(ifftshift(WF_image)));
    imwritestack(abs(WF_image_fft),[pathname 'SIM Result/Spectrum_Widefield-' filenameA(1:end-4) '.tif'])
end
clear temp WF_image

%% OPD estimation
Pindex = 1;
Phalf_index = 2;
if isOPLD == 1
    disp('With OPLD OTF reconstruction...')
    disp('Begin the registration of raw data A and B, please wait...')
    [SIM_raw6B,~]=Run_4ps_registration(SIM_rawA, SIM_rawB);
    disp('The registration of raw data A and B is completed.')
    clear SIM_rawA SIM_rawB
    % OTF resampling
    H=imreadstack(OTF_6B_Path);
    H_all = GenerateOTF_Fre_revised(H, Beads_pixesizeZ, pixelsizeZ, numz, max_xy);
    % OTF intensity normalisation
    if isNormalizOTF == 1
        H_all(:,:,1:numz) = H_all(:,:,1:numz)./max(max(max(H_all(:,:,1:numz)))); % m = 0
        H_all(:,:,1*numz+1:2*numz) = H_all(:,:,1*numz+1:2*numz)./max(max(max(H_all(:,:,1*numz+1:2*numz)))); % m = 2
        H_all(:,:,2*numz+1:3*numz) = H_all(:,:,2*numz+1:3*numz)./max(max(max(H_all(:,:,2*numz+1:3*numz)))); % m =1
    end
    H_zhongxin = H_all(:,:,Zmedian);                             % m = 0
    H_yiweiP = H_all(:,:,numz*(Pindex)+Zmedian);
    H_yiweiPhalf = H_all(:,:,numz*Phalf_index+Zmedian);          % m = 1
    pz=pixelsizeZ*numz*Pzrefmed/wavelengh;
    pz=round(pz);
    inteval=round(pz/2);
    pa=round(Zmedian-3/2*inteval); % OPD=a layer
    H_zhongxina = H_all(:,:,pa);                             %  m = 0
    H_yiweiPhalfa = H_all(:,:,numz*Phalf_index+pa);          %  m = 1
    clear H H_all

    %Parameters estimation
    zuobiaox = zeros(9,2);
    zuobiaoy = zeros(9,2);
    angle6 = zeros(nort,2);
    c6 = zeros(nort,2);
    angle = zeros(nort,1);
    disp('Fitting parameters, please wait...')
    for orti=1:nort
        [sep_im] = SIM_3D_seprt_PAZ(SIM_raw6B(:,:,1:endframe*nangle*nort),nangle,nort,spjg,regul,numz,orti); % [x,y,z,angle] Spectral separation of 5 phases in each direction
        % initial phase and contrast estimation
        spzhongxinAC = sep_im(:,:,Zmedian,3);
        spyiweiPAC = sep_im(:,:,Zmedian,Pindex);
        spyiweiPhalfAC = sep_im(:,:,Zmedian,Phalf_index);
        sep_im = sep_im./(abs(sep_im)+eps);
        spzhongxin = sep_im(:,:,Zmedian,3);
        spyiweiP = sep_im(:,:,Zmedian,Pindex);

        % OPD estimation
        spzhongxina = sep_im(:,:,pa,3);
        spyiweiPhalfa = sep_im(:,:,pa,Phalf_index);
        clear sep_im

        [tmpx,tmpy,tmp_angl,tmp_c6,~] = SIM_3D_p_CC_saveHistg(spzhongxin, spyiweiP, spzhongxinAC, spyiweiPAC, spyiweiPhalfAC,...
            H_zhongxin, H_yiweiP, H_yiweiPhalf,fc_para, pg, sizex, sizey, fanwei, orti, jindu); % Parameter estimation 2D correlation
        zuobiaox((orti-1)*3+1:(orti-1)*3+3,:) = tmpx((orti-1)*3+1:(orti-1)*3+3,:);
        zuobiaoy((orti-1)*3+1:(orti-1)*3+3,:) = tmpy((orti-1)*3+1:(orti-1)*3+3,:);
        angle6(orti,:) = tmp_angl;
        c6(orti,:) = tmp_c6;
        plong = ((zuobiaox((orti-1)*3+1:(orti-1)*3+3,:)-zuobiaox(1,1)).^2+ ...
            (zuobiaoy((orti-1)*3+1:(orti-1)*3+3,:)-zuobiaoy(1,1)).^2).^0.5;
        disp(['From ori' num2str(orti)]);
        disp(['plong: ' num2str(plong(2,:))]);

        [tmp_angl,~]=SIM_3D_p_CCa_saveHistg(angle6(orti,2),H_zhongxin,fc_para, sizex, sizey, orti, spzhongxina,spyiweiPhalfa,H_zhongxina,H_yiweiPhalfa,zuobiaox,zuobiaoy);
        angle(orti,:) = tmp_angl(1,2);
        opld(orti,:)=angle(orti,:)*wavelengh/Pzrefmed/2/pi;
        save([pathname 'SIM Result/' filenameA(1:end-4) '_parameter.mat'], 'zuobiaox', 'zuobiaoy', 'angle6', 'c6','opld');
    end

    OPLD_total=sum(opld,"all")./nort;
    disp('Contrast:')
    disp(c6)
    disp('Initial phase：')
    disp(angle6)
    disp('OPLD phase：')
    disp(angle)
%     disp('OPLD：')
%     disp(opld)
    disp('OPLD(nm)：')
    disp(OPLD_total)
    %% choose OPLD OTF: We provide two representative OTFs selected from the 4Pi-SIM OTF library.
    if 20 <= OPLD_total && OPLD_total < 30  % OPLD20
        H_real = imreadstack('OTF\4PSOTF_OPLD_20_real.tif');
        H_imag = imreadstack('OTF\4PSOTF_OPLD_20_imag.tif');
    elseif 180 <= OPLD_total && OPLD_total < 190  % OPLD180
        H_real = imreadstack('OTF\4PSOTF_OPLD_180_real.tif');
        H_imag = imreadstack('OTF\4PSOTF_OPLD_180_imag.tif');
    end

    HA=complex(H_real,H_imag);
    disp('Start to generate 3D OTF,please wait...')
    H_allA = GenerateOTF_Fre_revised(HA, Beads_pixesizeZ, pixelsizeZ, numz, max_xy);  % 3D OTF 疏采样, 三个相位的OTF 低频的OTF 移动p的OTF 移动p/2的OTF
    disp('3D OTF Generation is completed.')
    if isNormalizOTF == 1
        H_allA(:,:,1:numz) = H_allA(:,:,1:numz)./max(max(max(H_allA(:,:,1:numz)))); % m = 0
        H_allA(:,:,1*numz+1:2*numz) = H_allA(:,:,1*numz+1:2*numz)./max(max(max(H_allA(:,:,1*numz+1:2*numz)))); % m = 2
        H_allA(:,:,2*numz+1:3*numz) = H_allA(:,:,2*numz+1:3*numz)./max(max(max(H_allA(:,:,2*numz+1:3*numz)))); % m =1
    end
else  % dont estimate OPLD
    disp('Zero OPLD OTF reconstruction...')
    H_real = imreadstack('OTF\4PSOTF_OPLD_0_real.tif');
    H_imag = imreadstack('OTF\4PSOTF_OPLD_0_imag.tif');
    HA=complex(H_real,H_imag);
    % resampling
    disp('Start to generate 3D OTF,please wait...')
    H_allA = GenerateOTF_Fre_revised(HA, Beads_pixesizeZ, pixelsizeZ, numz, max_xy);
    disp('3D OTF Generation is completed.')
    if isNormalizOTF == 1
        H_allA(:,:,1:numz) = H_allA(:,:,1:numz)./max(max(max(H_allA(:,:,1:numz)))); % m = 0
        H_allA(:,:,1*numz+1:2*numz) = H_allA(:,:,1*numz+1:2*numz)./max(max(max(H_allA(:,:,1*numz+1:2*numz)))); % m = 2
        H_allA(:,:,2*numz+1:3*numz) = H_allA(:,:,2*numz+1:3*numz)./max(max(max(H_allA(:,:,2*numz+1:3*numz)))); % m =1
    end
    H_zhongxin = H_allA(:,:,Zmedian);                             %  m = 0
    H_yiweiP = H_allA(:,:,numz*(Pindex)+Zmedian);
    H_yiweiPhalf = H_allA(:,:,numz*Phalf_index+Zmedian);          %  m = 1
    clear HA
    % parameters estimation
    zuobiaox = zeros(9,2);
    zuobiaoy = zeros(9,2);
    angle6 = zeros(3,2);
    c6 = zeros(3,2);
    disp('Fitting parameters, please wait...')
    for orti=1:nort
        [sep_im] = SIM_3D_seprt_PAZ(SIM_rawA(:,:,1:endframe*nangle*nort),nangle,nort,spjg,regul,numz,orti);

        % initial phase and contrast
        spzhongxinAC = sep_im(:,:,Zmedian,3);             % m=0
        spyiweiPAC = sep_im(:,:,Zmedian,Pindex);          % m=2
        spyiweiPhalfAC = sep_im(:,:,Zmedian,Phalf_index); % m=1

        sep_im = sep_im./(abs(sep_im)+eps);
        spzhongxin = sep_im(:,:,Zmedian,3);
        spyiweiP = sep_im(:,:,Zmedian,Pindex);
        clear sep_im

        [tmpx,tmpy,tmp_angl,tmp_c6] = SIM_3D_p_CC_saveHistg(spzhongxin, spyiweiP, spzhongxinAC, spyiweiPAC, spyiweiPhalfAC,...
            H_zhongxin, H_yiweiP, H_yiweiPhalf,fc_para, pg, sizex, sizey, fanwei, orti, jindu); % Parameter estimation 2D correlation

        zuobiaox((orti-1)*3+1:(orti-1)*3+3,:) = tmpx((orti-1)*3+1:(orti-1)*3+3,:);
        zuobiaoy((orti-1)*3+1:(orti-1)*3+3,:) = tmpy((orti-1)*3+1:(orti-1)*3+3,:);
        angle6(orti,:) = tmp_angl;
        c6(orti,:) = tmp_c6;
    end
    save([pathname 'SIM Result/' filenameA(1:end-4) '_parameter.mat'], 'zuobiaox', 'zuobiaoy', 'angle6', 'c6');
    disp('Plong：')
    disp(((zuobiaox(2,1)-zuobiaox(1,1)).^2+(zuobiaoy(2,1)-zuobiaoy(1,1)).^2).^0.5)
    disp(((zuobiaox(5,1)-zuobiaox(1,1)).^2+(zuobiaoy(5,1)-zuobiaoy(1,1)).^2).^0.5)
    disp(((zuobiaox(8,1)-zuobiaox(1,1)).^2+(zuobiaoy(8,1)-zuobiaoy(1,1)).^2).^0.5)
    disp('Contrast:')
    disp(c6)
    disp('Initial phase：')
    disp(angle6)
end

clearvars -except pixelsize loadswitch H_allA H_allB pathname filenameA filenameB sizex sizey starframe numz nort nangle...
    bgflag bgname zuobiaox zuobiaoy angle6 c6 regul spjg weilac notch_swith notch_para ...
    pixelsize pixelsizeZ OTFmaskNA wavelengh Zpadding EmissionWavelength OTFmaskrefmed Pzrefmed notch_para_1 notch_para_2 ...
    isEqualization isSaveReconband I2M_para fc_delta M NA

% apodization parameters fc's value
fc_bhs = pixelsize*100/M*sizex*2*NA/EmissionWavelength+fc_delta;

% wiener deconvolution
% disp('Merging components,please wait...')
SIM_4PS_Wiener_deconv(H_allA,pathname,filenameA,sizex,sizey,numz,nort,nangle,...
    bgflag,bgname,zuobiaox,zuobiaoy,angle6,c6,regul,spjg,weilac,notch_swith,notch_para_1,notch_para_2,fc_bhs,...
    pixelsize,pixelsizeZ,wavelengh,Zpadding,isEqualization,I2M_para);

end
