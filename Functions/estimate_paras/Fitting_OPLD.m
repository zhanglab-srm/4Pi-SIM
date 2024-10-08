function [phase_opld] = Fitting_OPLD(FitBand, dataparams, Reconparams, orti)
% [phase_opld] = Fitting_OPLD(initial_phase,H0_para,ap_limit, sizex, sizey, orti, ...
%     band0_a,band1_a,H0_opld,H1_opld,coodinatex,coodinatey)
    % This function will return the estimated value of OPLD.

    hangs=2*orti-1;
    n_512 = max([dataparams.sizex,dataparams.sizey]);
    Pindex = 1;
    Phalf_index = 2;
    [otfx, otfy] = size(FitBand.H0_para);              
    [k_x_512, k_y_512] = meshgrid(-(n_512)/2:(n_512)/2-1, -(n_512)/2:(n_512)/2-1);
    k_r_512 = sqrt(k_x_512.^2+k_y_512.^2);     
    indi_512 =  k_r_512 > Reconparams.ap_limit ;             
    K_h = [otfx, otfy];
    N_h = 2*K_h;
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw = zeros(N_h);
    
    FitBand.H1_opld(indi_512)=0;
    FitBand.H0_opld(indi_512)=0;
    FitBand.H1_opld = abs(FitBand.H1_opld);
    FitBand.H0_opld = abs(FitBand.H0_opld);
    
    % pz_opld m=1 OTF
    H1_tmpopld = FitBand.H1_opld;              
    H1_tmpopld(H1_tmpopld~=0) = 1;
    hw(v_h{:}) = FitBand.H1_opld;             
    H1_opld = hw;
    hw(v_h{:}) = H1_tmpopld;
    H1_binary = hw;
    H1_opld = gpuArray(single(H1_opld));
    H1_binary = gpuArray(single(H1_binary));
    
    % pz_opld m=0 OTF
    H0_tmpopld = FitBand.H0_opld;                  
    H0_tmpopld(H0_tmpopld~=0) = 1;
    hw(v_h{:}) = FitBand.H0_opld;                  
    H0_opld = hw;
    hw(v_h{:}) = H0_tmpopld;
    H0_binary = hw;
    H0_opld = gpuArray(single(H0_opld));
    % H0_binary = gpuArray(single(H0_binary));
   
    % m=1 parameters
    maxx = Reconparams.coodinatex(hangs+1+floor((hangs-1)/2), Pindex);
    maxy = Reconparams.coodinatey(hangs+1+floor((hangs-1)/2), Pindex);
    [xx2,yy2] = meshgrid(-n_512:n_512-1, -n_512:n_512-1);
    xx2 = gpuArray(single(xx2));
    yy2 = gpuArray(single(yy2));
    maxx = n_512 + (maxx-n_512)/2;
    maxy = n_512 + (maxy-n_512)/2;
    kytest = 2*pi*(maxx-n_512)/(2*n_512);
    kxtest = 2*pi*(maxy-n_512)/(2*n_512);
    Irtest = exp(1i*(kxtest*xx2+kytest*yy2));

    hw(v_h{:}) = FitBand.band0_a;    
    band0_a = hw;
    hw(v_h{:}) = FitBand.band1_a;   
    band1_a = hw;

    % estimate opld
    phasea = exp(-1i*Reconparams.initial_phase(orti,2));    % remove phase a from spectrum (see Supplementary Figure 6i,j and Methods)
    moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_binary))).*Irtest)));
    moveMask = abs(moveMask);
    moveMask(abs(moveMask)>0.9) = 1;
    moveMask(abs(moveMask)~=1) = 0;
    moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_opld))).*Irtest))));
    moveOTF = moveOTF.*moveMask;
    temp_band = gather(fftshift(fft2(ifft2(ifftshift(band1_a)).*Irtest)).*moveMask.*H0_opld.*phasea);
    
    [phase_opld(Phalf_index),~] = Fitting_histogram(gather(band0_a),gather(moveOTF),temp_band,gather(H0_opld));

end