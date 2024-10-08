function [Irtest, moveMask, maxx, maxy, moveOTF, H0_OTF, H1_OTF, H1_binary,hw,v_h,xx2,yy2] = Fitting_plong(ap_limit,p_theory,p_limit,p_precision, ...
    sizex,sizey,Pindex,max_xy,FitBand)
    % This function will estimate the wave vector for structured illumination
    % The integer-pixel coordinates of the illumination wave vector k0 are estimated by cross-correlating the shifted order-2 and the order-0 SIM bands.
    % Subpixel precision in locating k0 is achieved by utilizing gradient descent methodology.

    %% initialization
    [otfx, otfy] = size(FitBand.H0_para);                  
    [k_x_512, k_y_512] = meshgrid(-(max_xy)/2:(max_xy)/2-1, -(max_xy)/2:(max_xy)/2-1);
    k_r_512 = sqrt(k_x_512.^2 + k_y_512.^2);    
    indi_512 =  k_r_512 > ap_limit ;              
    K_h = [otfx, otfy];
    N_h = 2*K_h;
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw = zeros(N_h);
    FitBand.H0_para(indi_512) = 0;                  
    FitBand.H2_para(indi_512) = 0;
    FitBand.H1_para(indi_512) = 0;
    FitBand.H0_para = abs(FitBand.H0_para);
    FitBand.H2_para = abs(FitBand.H2_para);
    FitBand.H1_para = abs(FitBand.H1_para);
    % m=0 OTF
    H0_temp = FitBand.H0_para;                  
    H0_temp(H0_temp~=0) = 1;
    hw(v_h{:}) = FitBand.H0_para;                   
    H0_OTF = hw;
    hw(v_h{:}) = H0_temp;
    H0_binary = hw;
    H0_OTF = gpuArray(single(H0_OTF));
    H0_binary = gpuArray(single(H0_binary));
    % m=2 OTF
    H2_temp = FitBand.H2_para;                        
    H2_temp(H2_temp~=0) = 1;
    hw(v_h{:}) = FitBand.H2_para;                    
    H2_OTF = hw;
    hw(v_h{:}) = H2_temp;
    H2_binary = hw;
    H2_OTF = gpuArray(single(H2_OTF));
    H2_binary = gpuArray(single(H2_binary));
    % m=1 OTF
    H1_temp = FitBand.H1_para;             
    H1_temp(H1_temp~=0) = 1;
    hw(v_h{:}) = FitBand.H1_para;                 
    H1_OTF = hw;
    hw(v_h{:}) = H1_temp;
    H1_binary = hw;
    H1_OTF = gpuArray(single(H1_OTF));
    H1_binary = gpuArray(single(H1_binary));
    clear H0_temp H1_temp H2_temp 
    % H0_para H1_para H2_para
    
    %% calculation
    stepx = 1;
    stepy = 1;

    % use cross-correlation to find integer plong coodinates
    Hnorm = ones(sizex, sizey);
    Normcoeff = dft(Hnorm,Hnorm);               
    sp_tmp = conj(FitBand.band2_p);
    sp_tmp = sp_tmp(end:-1:1,end:-1:1);      
    H_CC = dft(FitBand.band0_p,sp_tmp);         
    plong_int = abs(H_CC./(Normcoeff+eps));
    [k_x_1024, k_y_1024] = meshgrid(-(2*max_xy)/2+1:(2*max_xy)/2-1, -(2*max_xy)/2+1:(2*max_xy)/2-1);
    k_r_1024 = sqrt(k_x_1024.^2+k_y_1024.^2);
    indi = ((k_r_1024 < (floor(p_theory/Pindex)-p_limit))|(k_r_1024 > (floor(p_theory/Pindex)+p_limit)));  
    plong_int(indi) = 0;
    [xx,yy] = ind2sub(size(plong_int),find(plong_int == max(plong_int(:)))); 
  
    % use gradient descent to find sub-pixel plong coodinates
    maxx = xx;
    maxy = yy;
    ky = 2*pi*(maxx-max_xy)/(2*max_xy);
    kx = 2*pi*(maxy-max_xy)/(2*max_xy);
    [xx2,yy2] = meshgrid(-max_xy:max_xy-1, -max_xy:max_xy-1);  
    xx2 = gpuArray(single(xx2));
    yy2 = gpuArray(single(yy2));
    hw(v_h{:}) = FitBand.band2_p;    
    band2_p = hw;
    band2_p = gpuArray(single(band2_p));
    hw(v_h{:}) = FitBand.band0_p;    
    band0_p = hw;
    band0_p = gpuArray(single(band0_p));

    Irtest = exp(1i*(kx*xx2+ky*yy2));
    moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_binary))).*Irtest))); 
    moveMask = abs(moveMask);
    moveMask(abs(moveMask)>0.9) = 1; 
    moveMask(abs(moveMask)~=1) = 0;
    moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_OTF))).*Irtest))));  
    moveOTF = moveOTF.*moveMask; 
    overlapRegion = (conj(fftshift(fft2(ifft2(ifftshift(band2_p)).*Irtest)).*moveMask.*H0_OTF)).*(band0_p.*moveOTF);
    overlapRegion = abs(sum(overlapRegion(:)));
    normRegion = H0_binary.*moveMask; 
    normRegion = sum(normRegion(:));
    overlapRegion = overlapRegion./normRegion;
    he = overlapRegion;
    % determine the direction of the gradient
    maxx_tmp1 = maxx-10^-5;
    maxx_tmp2 = maxx+10^-5;
    maxy_tmp1 = maxy-10^-5;
    maxy_tmp2 = maxy+10^-5;
    ky_tmp1 = 2*pi*(maxx_tmp1-max_xy)/(2*max_xy);
    kx_tmp1 = 2*pi*(maxy_tmp1-max_xy)/(2*max_xy);
    ky_tmp2 = 2*pi*(maxx_tmp2-max_xy)/(2*max_xy);
    kx_tmp2 = 2*pi*(maxy_tmp2-max_xy)/(2*max_xy);
    for ii=1:1:4  % find the optimum in different directions
        switch ii
            case 1
                kxtest=kx_tmp1;
                kytest=ky;
            case 2
                kxtest=kx_tmp2;
                kytest=ky;
            case 3
                kxtest=kx;
                kytest=ky_tmp1;
            case 4
                kxtest=kx;
                kytest=ky_tmp2;
        end
        Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
        moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_binary))).*Irtest)));
        moveMask = abs(moveMask);
        moveMask(abs(moveMask)>0.9) = 1;
        moveMask(abs(moveMask)~=1) = 0;
        moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_OTF))).*Irtest))));
        moveOTF = moveOTF.*moveMask;
        overlapRegion = (conj(fftshift(fft2(ifft2(ifftshift(band2_p)).*Irtest)).*moveMask.*H0_OTF)).*(band0_p.*moveOTF);
        overlapRegion = abs(sum(overlapRegion(:)));
        normRegion = H0_binary.*moveMask;
        normRegion = sum(normRegion(:));
        overlapRegion = overlapRegion./normRegion;
        test(ii) = overlapRegion;
    end
    if((test(1)>test(2)))
        flag_maxy=-1;
    elseif((test(1)<test(2)))
        flag_maxy=+1;
    else
        flag_maxy=+1;
    end
    if((test(3)>test(4)))
        flag_maxx=-1;
    elseif((test(3)<test(4)))
        flag_maxx=+1;
    else
        flag_maxx=+1;
    end
    while ((stepx > p_precision) || (stepy > p_precision))
        %%%%%%%%%%%%%%%%%%%%%%% maxxorientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        maxx_tmp1 = maxx-10^-5;
        maxx_tmp2 = maxx+10^-5;
        ky_tmp1 = 2*pi*(maxx_tmp1-max_xy)/(2*max_xy);
        ky_tmp2 = 2*pi*(maxx_tmp2-max_xy)/(2*max_xy);
        for ii=3:1:4   
            switch ii
                case 3
                    kxtest = 2*pi*(maxy-max_xy)/(2*max_xy);
                    kytest = ky_tmp1;
                case 4
                    kxtest = 2*pi*(maxy-max_xy)/(2*max_xy);
                    kytest = ky_tmp2;
            end
            Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
            moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_binary))).*Irtest))); 
            moveMask = abs(moveMask);
            moveMask(abs(moveMask)>0.9) = 1;
            moveMask(abs(moveMask)~=1) = 0;
            moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_OTF))).*Irtest)))); 
            moveOTF = moveOTF.*moveMask;  
            overlapRegion = (conj(fftshift(fft2(ifft2(ifftshift(band2_p)).*Irtest)).*moveMask.*H0_OTF)).*(band0_p.*moveOTF); 
            overlapRegion = abs(sum(overlapRegion(:)));
            normRegion = H0_binary.*moveMask;
            normRegion = sum(normRegion(:));
            overlapRegion = overlapRegion./normRegion;
            test(ii) = overlapRegion;
        end
        if((test(3)>test(4)))
            flag_maxx=-1;
        elseif((test(3)<test(4)))
            flag_maxx=+1;
        else
            flag_maxx=-1* flag_maxx;
        end
        while(stepx>(p_precision))
            maxx_tmp = maxx+flag_maxx*stepx;
            kytest = 2*pi*(maxx_tmp-max_xy)/(2*max_xy);
            kxtest = 2*pi*(maxy-max_xy)/(2*max_xy);
            Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
            moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_binary))).*Irtest)));
            moveMask = abs(moveMask);
            moveMask(abs(moveMask)>0.9)=1;
            moveMask(abs(moveMask)~=1)=0;
            moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_OTF))).*Irtest))));
            moveOTF = moveOTF.*moveMask;
            overlapRegion = (conj(fftshift(fft2(ifft2(ifftshift(band2_p)).*Irtest)).*moveMask.*H0_OTF)).*(band0_p.*moveOTF);
            overlapRegion = abs(sum(overlapRegion(:)));
            normRegion = H0_binary.*moveMask;
            normRegion = sum(normRegion(:));
            overlapRegion = overlapRegion./normRegion;
            he_tmp = overlapRegion;
            if(he_tmp <= he)
                stepx = 0.5*stepx;
            elseif(he_tmp > he)
                he = he_tmp;      
                maxx = maxx_tmp;  
                break;
            end
        end    % end of while x
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% maxyorientation %%%%%%%%%%%%%%%%%%%%
        maxy_tmp1 = maxy-10^-5;
        maxy_tmp2 = maxy+10^-5;
        kx_tmp1 = 2*pi*(maxy_tmp1-max_xy)/(2*max_xy);
        kx_tmp2 = 2*pi*(maxy_tmp2-max_xy)/(2*max_xy);
        for ii=1:1:2 
            switch ii
                case 1
                    kxtest = kx_tmp1;
                    kytest = 2*pi*(maxx-max_xy)/(2*max_xy);
                case 2
                    kxtest = kx_tmp2;
                    kytest = 2*pi*(maxx-max_xy)/(2*max_xy);
            end
            Irtest=exp(1i*(kxtest*xx2+kytest*yy2));
            moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_binary))).*Irtest)));
            moveMask = abs(moveMask);
            moveMask(abs(moveMask)>0.9)=1;
            moveMask(abs(moveMask)~=1)=0;
            moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_OTF))).*Irtest))));
            moveOTF= moveOTF.*moveMask;
            overlapRegion = (conj(fftshift(fft2(ifft2(ifftshift(band2_p)).*Irtest)).*moveMask.*H0_OTF)).*(band0_p.*moveOTF);
            overlapRegion = abs(sum(overlapRegion(:)));
            normRegion = H0_binary.*moveMask;
            normRegion = sum(normRegion(:));
            overlapRegion = overlapRegion./normRegion;
            test(ii) = overlapRegion;
        end
        if((test(1)>test(2)))
            flag_maxy=-1;
        elseif((test(1)<test(2)))
            flag_maxy=+1;
        else
            flag_maxy=-1*flag_maxy;
        end
        while(stepy > (p_precision))
            maxy_tmp = maxy+flag_maxy*stepy;
            kytest = 2*pi*(maxx-max_xy)/(2*max_xy);
            kxtest = 2*pi*(maxy_tmp-max_xy)/(2*max_xy);
            Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
            moveMask = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_binary))).*Irtest)));
            moveMask = abs(moveMask);
            moveMask(abs(moveMask)>0.9) = 1;
            moveMask(abs(moveMask)~=1) = 0;
            moveOTF = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H2_OTF))).*Irtest))));
            moveOTF = moveOTF.*moveMask;
            overlapRegion = (conj(fftshift(fft2(ifft2(ifftshift(band2_p)).*Irtest)).*moveMask.*H0_OTF)).*(band0_p.*moveOTF);
            overlapRegion = abs(sum(overlapRegion(:)));
            normRegion = H0_binary.*moveMask;
            normRegion = sum(normRegion(:));
            overlapRegion = overlapRegion./normRegion;
            he_tmp = overlapRegion;
            if(he_tmp <= he)
                stepy = 0.5*stepy;
            elseif(he_tmp > he)
                he = he_tmp;    
                maxy = maxy_tmp;
                break;
            end
        end       % end of while y
    end
    
end