function Wiener_deconv_4PS(dataparams, OTFparams, Reconparams)
    % This function performs Wiener filtering to SIM data.
    disp('Performing Wiener filtering...');
    % prepare OTF
    OTFparams.H_Recon = single(OTFparams.H_Recon);
    H0_OTF = OTFparams.H_Recon(:,:,1:dataparams.numz);            
    H1_OTF = OTFparams.H_Recon(:,:,2*dataparams.numz+1:dataparams.numz*3); 
    H2_OTF = OTFparams.H_Recon(:,:,1*dataparams.numz+1:dataparams.numz*2);
    clear Reconparams.H_Recon
    OTF_xy = size(H0_OTF,1);
    OTF_zpad = dataparams.numz+2*Reconparams.Zpadding;
    H0_OTF = imresize3(H0_OTF, [OTF_xy, OTF_xy, OTF_zpad], 'linear');
    H1_OTF = imresize3(H1_OTF, [OTF_xy, OTF_xy, OTF_zpad], 'linear');
    H2_OTF = imresize3(H2_OTF, [OTF_xy, OTF_xy, OTF_zpad], 'linear');
    all_name = [dataparams.pathname,dataparams.filenameA];
    
    % background noise
    if Reconparams.bgflag == 1
        bg = single(imreadstack(Reconparams.bgname));
        [bgsizex, bgsizey] = size(bg);
        if bgsizex == dataparams.sizex && bgsizey == dataparams.sizey
        else
            bg = ones(dataparams.sizex, dataparams.sizey, 'single').*100;
        end
    elseif Reconparams.bgflag == 0
        bg = zeros(dataparams.sizex,dataparams.sizey, 'single');
    end
    
    % rescale coordinate
    coodinatex = Reconparams.coodinatex*(OTF_xy/Reconparams.coodinatex(1,1));
    coodinatey = Reconparams.coodinatey*(OTF_xy/Reconparams.coodinatey(1,1));
    
    % load parameters
    contrast = Reconparams.contrast;
    modulate_depth = [(contrast(1,1));(contrast(1,2));1;(contrast(1,2));(contrast(1,1));
        (contrast(2,1));(contrast(2,2));1;(contrast(2,2));(contrast(2,1));
        (contrast(3,1));(contrast(3,2));1;(contrast(3,2));(contrast(3,1));];
    ap_range = ceil(Reconparams.ap_range*(OTF_xy/512));  % maximum spatial frequency
    plong = sum(((coodinatex-coodinatex(1,1)).^2+(coodinatey-coodinatey(1,1)).^2).^0.5, 1)./dataparams.nort./2;
    movebands = zeros((2*OTF_xy), (2*OTF_xy), OTF_zpad, 'single'); 
    tempnum1 = dataparams.nort*dataparams.nangle;
    tempnum2 = dataparams.nort*dataparams.nangle*dataparams.numz;
    Irtest = zeros((2*OTF_xy), (2*OTF_xy), tempnum1, 'single');      
    [xx2, yy2] = meshgrid(-OTF_xy:OTF_xy-1, -OTF_xy:OTF_xy-1);                    
    [k_x, k_y] = meshgrid(-(OTF_xy)/2:(OTF_xy)/2-1, -(OTF_xy)/2:(OTF_xy)/2-1);
    k_r = sqrt(k_x.^2+k_y.^2);
    clear k_x k_y
    
    % transform coordinates
    [tmp_coodinatex,tmp_coodinatey] = Transform_coodinate(Reconparams.coodinatex,Reconparams.coodinatey);  
    
    % prepare OTF mask
    Mask_4PS = imreadstack('OTF/Mask_4PS.tif');
    Mask_4PS = Resample_OTF(Mask_4PS, dataparams.pixelsizeZ, dataparams.pixelsizeZ, dataparams.numz, OTF_xy);
    otfmask0 = Mask_4PS(:,:,1:size(Mask_4PS,3)/3);
    otfmask2 = Mask_4PS(:,:,(size(Mask_4PS,3)/3+1):2*size(Mask_4PS,3)/3);
    otfmask1 = Mask_4PS(:,:,(2*size(Mask_4PS,3)/3+1):size(Mask_4PS,3));
    clear Mask_4PS
    otfmask0 = imresize3(single(otfmask0), [OTF_xy, OTF_xy, OTF_zpad], 'linear');
    otfmask0(otfmask0<0.1) = 0;
    otfmask0(otfmask0~=0) = 1;
    otfmask1 = imresize3(single(otfmask1), [OTF_xy, OTF_xy, OTF_zpad], 'linear');
    otfmask1(otfmask1<0.1) = 0;
    otfmask1(otfmask1~=0) = 1;
    otfmask2 = imresize3(single(otfmask2), [OTF_xy, OTF_xy, OTF_zpad], 'linear');
    otfmask2(otfmask2<0.1) = 0;
    otfmask2(otfmask2~=0) = 1;
    H0OTF_Pad = H0_OTF.*otfmask0;       % beyond the cut-off frequency is 0
    H2OTF_Pad = H2_OTF.*otfmask2;       
    H1OTF_Pad = H1_OTF.*otfmask1;      
    
    % padding
    K_h = size(H0OTF_Pad);
    if numel(K_h) == 3
        N_h = 2*K_h-[0,0,K_h(1,3)];
    else   
        N_h = 2*K_h-[0,0];
    end
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw = zeros(N_h, 'single');
    hw(v_h{:}) = H0OTF_Pad;
    H0OTF_Pad = single(hw);
    hw(v_h{:}) = H2OTF_Pad;
    H2OTF_Pad = single(hw);
    hw(v_h{:}) = H1OTF_Pad;
    H1OTF_Pad = single(hw);
    clear hw N_h K_h L_h v_h
    
    % prepare notch filter
    if Reconparams.notch_switch == 1
        notch_filter = 1-exp(-Reconparams.notch_para_1*(abs(k_r).^Reconparams.notch_para_2));
        notch_filter = single(repmat(notch_filter,[1,1,dataparams.nangle]));
    elseif Reconparams.notch_switch == 0
        notch_filter = ones(size(k_r), 'single');
        notch_filter = repmat(notch_filter,[1,1,dataparams.nangle]);
    end
    notch_filter = single(notch_filter);
    
    % prepare shift matrix
    for ii=1:tempnum1
        kytest = 2*pi*(tmp_coodinatex(ii)-OTF_xy)/(2*OTF_xy);
        kxtest = 2*pi*(tmp_coodinatey(ii)-OTF_xy)/(2*OTF_xy);
        Irtest(:,:,ii) = exp(1i*(kxtest*xx2+kytest*yy2)); 
    end
    
    % prepare apodization function
    [apod_func] = apodization(OTF_xy, dataparams.numz, plong, ap_range, ...
        dataparams.pixelsizeZ, dataparams.Exwave); 
    apod_func = single(imresize3(apod_func, [2*OTF_xy, 2*OTF_xy, OTF_zpad], 'linear'));
    
    % prepare Wiener denominator
    disp('  Calculating Wiener denominator...');
    WienerDenom = zeros((2*OTF_xy), (2*OTF_xy), OTF_zpad,  'single');
    for ii = 1:tempnum1
        if mod(ii,dataparams.nangle)==ceil(dataparams.nangle/2)-1 || mod(ii,dataparams.nangle)==ceil(dataparams.nangle/2)+1 
            movebands = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(modulate_depth(ii).* ...
                H1OTF_Pad))).*Irtest(:,:,ii)) )); % OTF  m=1
            OTFmask=otfmask1;
        elseif mod(ii,dataparams.nangle)==ceil(dataparams.nangle/2)
            movebands = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(modulate_depth(ii).* ...
                H0OTF_Pad))).*Irtest(:,:,ii)) )); % OTF  m=0
            OTFmask=otfmask0;
        else
            movebands = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(modulate_depth(ii).* ...
                H2OTF_Pad))).*Irtest(:,:,ii)) )); % OTF  m=2
            OTFmask=otfmask2;
        end
    
        OTFmask = imresize3(single(OTFmask), [2*OTF_xy, 2*OTF_xy, OTF_zpad], 'linear');
        OTFmask = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(OTFmask))).*Irtest(:,:,ii)) ));
        OTFmask(OTFmask<0.1) = 0;
        OTFmask(OTFmask~=0) = 1;
    
        movebands = abs(movebands);
        movebands = movebands.*OTFmask; 
        WienerDenom = WienerDenom + movebands.^2;
    end
    clear movebands OTFmask otfmask1 otfmask0 otfmask2
    Denom_Shift = zeros(OTF_xy, OTF_xy, OTF_zpad, tempnum1, 'single');
    apod_Shift = zeros(OTF_xy, OTF_xy, OTF_zpad, tempnum1, 'single');
    for ii= 1:tempnum1
        Denom_Shift(:,:,:,ii) = Crop(fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift( WienerDenom ))).* ...
            conj(Irtest(:,:,ii))))),H0_OTF);
        apod_Shift(:,:,:,ii) = Crop(fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift( apod_func ))).* ...
            conj(Irtest(:,:,ii))))),H0_OTF);
    end
    WienerDenom_Move = abs(apod_Shift)./(abs(Denom_Shift) + 0.005*(Reconparams.wiener_coeff)^2);
    clear WienerDenom Denom_Shift apod_Shift sigma x y
    
    % prepare separation phase matrix
    [phase_matrix] = Make_phase_matrix(Reconparams.initial_phase,dataparams.spacing,dataparams.regul);
    
    % make mask and window function
    padsize = 0;
    x = 1:(dataparams.sizey+2*padsize);
    y = (1:(dataparams.sizex+2*padsize))';
    sigma = 0.25;
    mask = single(repmat(sigmoid(sigma*(x-padsize)) - sigmoid(sigma*(x-dataparams.sizey-padsize-1)), ...
        dataparams.sizex+2*padsize, 1).* repmat(sigmoid(sigma*(y-padsize)) - sigmoid(sigma* ...
        (y-dataparams.sizex-padsize-1)), 1, dataparams.sizey+2*padsize));

    K_h2 = [OTF_xy, OTF_xy, OTF_zpad];
    N_h2 = [2*OTF_xy, 2*OTF_xy, OTF_zpad];
    L_h2 = ceil((N_h2-K_h2) / 2);
    v_h2 = colonvec(L_h2+1, L_h2+K_h2);
    K_h1 = [dataparams.sizex,dataparams.sizey,OTF_zpad];
    N_h1 = [OTF_xy,OTF_xy,OTF_zpad];
    L_h1 = ceil((N_h1-K_h1) / 2);
    v_h1 = colonvec(L_h1+1, L_h1+K_h1);
    hw1 = zeros(N_h1, 'single');
    wd = Window(OTF_zpad, dataparams.sizex, dataparams.sizey, 0.25, 3);
    hw1(v_h1{:}) = wd;
    wd = hw1;
        
    clearvars -except Reconparams dataparams OTFparams all_name bg notch_filter mask WienerDenom_Move Irtest wd ...
        apod_func OTF_xy OTF_zpad phase_matrix modulate_depth v_h2 N_h2 H0_OTF H1_OTF H2_OTF tempnum1 tempnum2
    
    sep_bands = zeros(OTF_xy, OTF_xy, OTF_zpad, tempnum1, 'single');
    reconimage = zeros(2*OTF_xy, 2*OTF_xy, OTF_zpad, 'single');
    num = 0;
    num = num + 1;
    if ( 1+(num-1)*tempnum2 + tempnum2-1 <= tempnum2 )
        SIM_raw = single(myimreadstack_TIRF(all_name, 1+(num-1)*tempnum2, tempnum2, ...
            dataparams.sizex, dataparams.sizey));
    else
        error('Length overflows');
    end
    SIM_raw = SIM_raw - bg;
    SIM_raw(SIM_raw<0) = 0;
    
    % balance image intensity
    restack = zeros(dataparams.sizex,dataparams.sizey,dataparams.numz,tempnum1);
    weights = zeros(1,dataparams.nangle*dataparams.nort);
    for ii = 1:dataparams.nort
        for jj=1:dataparams.nangle
            restack(:,:,1:dataparams.numz,(ii-1)*dataparams.nangle+jj) = SIM_raw(:,:,(ii-1)* ...
                dataparams.nangle+jj:tempnum1:tempnum2);
            weights(:,(ii-1)*dataparams.nangle+jj) = mean(mean(mean(restack(:,:,:, ...
                (ii-1)*dataparams.nangle+jj))));
        end
    end
    clear restack 
    weights = mean(weights)./weights;

    % separate bands
    disp('  Separating raw data...');
    SIM_raw = SIM_raw.*(mask.^3);
    Zpadding = Reconparams.Zpadding;
    for orti=1:dataparams.nort
        phase_im = zeros(size(SIM_raw,1),size(SIM_raw,2),OTF_zpad,dataparams.nangle,'single');
        for anglei=1:dataparams.nangle 
            phase_im(:,:,Zpadding+1:end-Zpadding,anglei) = SIM_raw(:,:,anglei+(orti-1)* ...
                dataparams.nangle:tempnum1:tempnum2).*weights(:,anglei+(orti-1)*dataparams.nangle);
            if Zpadding~=0
                phase_im(:,:,1:Zpadding,anglei) = repmat(SIM_raw(:,:,anglei+(orti-1)* ...
                    dataparams.nangle),[1,1,Zpadding]);
                phase_im(:,:,end-Zpadding+1:end,anglei) = repmat(SIM_raw(:,:,tempnum2-(tempnum1-(anglei+(orti-1)* ...
                    dataparams.nangle))),[1,1,Zpadding]);
            end
            phase_im(:,:,:,anglei) = phase_im(:,:,:,anglei).*wd;
        end
        ztoxy = zeros(size(phase_im,4),size(phase_im,1)*size(phase_im,2)*size(phase_im,3),'single');
        for testi=1:size(phase_im,4)
            temp = (phase_im(:,:,:,testi));
            ztoxy(testi,:) = temp(:);
        end
        clear temp
        ztoxy = phase_matrix(:,:,orti) * ztoxy;
        for testi=1:size(phase_im,4)  
            sep_bands(:,:,:,(orti-1)*dataparams.nangle + testi) = reshape(ztoxy(testi,:), ...
                [size(phase_im,1),size(phase_im,2),size(phase_im,3)]); 
            sep_bands(:,:,:,(orti-1)*dataparams.nangle + testi) = fftshift(fftn(sep_bands(:,:,:,(orti-1)* ...
                dataparams.nangle + testi)));
            sep_bands(:,:,:,(orti-1)*dataparams.nangle + testi) = sep_bands(:,:,:,(orti-1)* ...
                dataparams.nangle + testi).* notch_filter(:,:,testi);
        end
    end
    clear ztoxy phase_im SIM_raw sep_im
    
    % recombine filtered bands
    disp('  Recombining filtered bands...');
    hw2 = zeros(N_h2,'single');
    for t = 1:tempnum1
        if mod(t,dataparams.nangle) == ceil(dataparams.nangle/2)-1 || mod(t,dataparams.nangle)==ceil(dataparams.nangle/2)+1
            sep_bands(:,:,:,t) = modulate_depth(t).* sep_bands(:,:,:,t).* conj(H1_OTF).* WienerDenom_Move(:,:,:,t);
        elseif mod(t,dataparams.nangle) == ceil(dataparams.nangle/2)
            sep_bands(:,:,:,t) = modulate_depth(t).* sep_bands(:,:,:,t).* conj(H0_OTF).* WienerDenom_Move(:,:,:,t);
        else
            sep_bands(:,:,:,t) = modulate_depth(t).* sep_bands(:,:,:,t).* conj(H2_OTF).* WienerDenom_Move(:,:,:,t);
        end
        hw2(v_h2{:}) = sep_bands(:,:,:,t);
        reconimage = reconimage + ifftn(ifftshift(hw2)).*Irtest(:, :, t);
    end
    clear sep_bands hw2
    reconimage = real(reconimage);
    reconimage = reconimage((end/2)+1-(dataparams.sizex):(end/2)+(dataparams.sizex),(end/2)+1- ...
        (dataparams.sizey):(end/2)+(dataparams.sizey),:);
    reconimage(reconimage<0) = 0;
    
    % save data
    if num == 1
        imwritestack(abs(fftshift(fftn(reconimage(:,:,Zpadding+1:end-Zpadding)))), [dataparams.pathname ...
            'SIM Result/SpectrumWiener' dataparams.filenameA(1:end-4) '.tif']);
    else
        imwritestacka(abs(fftshift(fftn(reconimage(:,:,Zpadding+1:end-Zpadding)))), [dataparams.pathname ...
            'SIM Result/SpectrumWiener' dataparams.filenameA(1:end-4) '.tif']);
    end
    if num == 1
        imwritestack(reconimage(:,:,Zpadding+1:end-Zpadding), [dataparams.pathname 'SIM Result/Wiener' ...
            dataparams.filenameA(1:end-4) '.tif']);
    else
        imwritestacka(reconimage(:,:,Zpadding+1:end-Zpadding), [dataparams.pathname 'SIM Result/Wiener' ...
            dataparams.filenameA(1:end-4) '.tif']);
    end
    clear reconimage
    disp('The 4Pi-SIM image reconstruction is completed.');
    disp('Please see SIM Result file for the reconstruction images.');
end