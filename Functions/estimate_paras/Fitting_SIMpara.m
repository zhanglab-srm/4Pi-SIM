function [coodinatetmpx,coodinatetmpy,tmpphase,tmpc] = Fitting_SIMpara(FitBand, dataparams, Reconparams, orti)
% [coodinatetmpx,coodinatetmpy,tmpphase,tmpc] = Fitting_SIMpara(band0_p, band2_p, band0_para, band2_para, band1_para,...
%             H0_para, H2_para, H1_para,ap_limit, p_theory, p_limit, p_precision, sizex, sizey, orti)
    % This function will return the estimated value of initial phase and contrast.

    hangs = 2*orti-1;
    max_xy = max([dataparams.sizex,dataparams.sizey]);
    n = 512;
    p_theory = Reconparams.p_theory*(max_xy/n);
    ap_limit = Reconparams.ap_limit*(max_xy/n);
    Pindex = 1;
    Phalf_index = 2;
    coodinatetmpx(1:9,1:2) = max_xy;              
    coodinatetmpy(1:9,1:2) = max_xy; 

    [Irtest, moveMask, maxx, maxy, moveOTF, H0_OTF, H1_OTF, H1_binary,hw,v_h,xx2,yy2] = Fitting_plong(ap_limit,p_theory, ...
        Reconparams.p_limit,Reconparams.p_precision, dataparams.sizex,dataparams.sizey,Pindex,max_xy,FitBand);

    hw(v_h{:}) = FitBand.band2_para;       
    band2_para = hw;
    hw(v_h{:}) = FitBand.band0_para;   
    band0_para = hw;
    hw(v_h{:}) = FitBand.band1_para;   
    band1_para = hw;
    
    % m=2 fitting parameters
    temp_band = gather(fftshift(fft2(ifft2(ifftshift(band2_para)).*Irtest)).*moveMask.*H0_OTF);
    coodinatetmpx(hangs+1+floor((hangs-1)/2), Pindex) = maxx;           
    coodinatetmpy(hangs+1+floor((hangs-1)/2), Pindex) = maxy;
    coodinatetmpx(hangs+2+floor((hangs-1)/2), Pindex) = 2*max_xy-maxx;   
    coodinatetmpy(hangs+2+floor((hangs-1)/2), Pindex) = 2*max_xy-maxy;
    [tmpphase(Pindex), tmpc(Pindex)] = Fitting_histogram(gather(band0_para),gather(moveOTF),temp_band,gather(H0_OTF));
    
    % m=1 fitting parameters
    maxx = max_xy + (maxx-max_xy)/2;
    maxy = max_xy + (maxy-max_xy)/2;
    kytest = 2*pi*(maxx-max_xy)/(2*max_xy);
    kxtest = 2*pi*(maxy-max_xy)/(2*max_xy);
    Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
    moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_binary))).*Irtest)));
    moveMask = abs(moveMask);
    moveMask(abs(moveMask)>0.9) = 1;
    moveMask(abs(moveMask)~=1) = 0;
    moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_OTF))).*Irtest))));
    moveOTF = moveOTF.*moveMask;
    temp_band = gather(fftshift(fft2(ifft2(ifftshift(band1_para)).*Irtest)).*moveMask.*H0_OTF);
    coodinatetmpx(hangs+1+floor((hangs-1)/2),Phalf_index) = maxx;          
    coodinatetmpy(hangs+1+floor((hangs-1)/2),Phalf_index) = maxy;
    coodinatetmpx(hangs+2+floor((hangs-1)/2),Phalf_index) = 2*max_xy-maxx;   
    coodinatetmpy(hangs+2+floor((hangs-1)/2),Phalf_index) = 2*max_xy-maxy;
    
    [tmpphase(Phalf_index),tmpc(Phalf_index)] = Fitting_histogram(gather(band0_para),gather(moveOTF),temp_band,gather(H0_OTF));
    
end