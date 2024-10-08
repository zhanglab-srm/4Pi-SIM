function  [sep_im] = Separate_SIMbands(SIM_raw,nangle,nort,spacing,regul,numz,orti)
    % This function will separate SIM spectrum bands.

    n = max([size(SIM_raw,1),size(SIM_raw,2)]);
    K_h = [size(SIM_raw,1),size(SIM_raw,2)];
    N_h = [n, n];
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw = zeros(N_h);
    temp1 = zeros(n, n, size(SIM_raw,3));
    for ii=1:1:size(SIM_raw,3)  
        hw(v_h{:}) = SIM_raw(:,:,ii);
        temp1(:,:,ii) = hw;
    end
    SIM_raw = temp1;
    clear temp1 hw;
    phase_im = zeros(size(SIM_raw,1), size(SIM_raw,2), numz, nangle);
    for anglei=1:nangle  
        phase_im(:,:,:,anglei) = SIM_raw(:,:,anglei+(orti-1)*nangle:nangle*nort:numz*nangle*nort);
    end
    clear SIM_raw
    phase_matrix= [1 1 1 1 1; % separation phase matrix
        exp(2*1i*regul*(spacing(1)/sum(spacing)))  exp(1i*regul*(spacing(1)/sum(spacing))) 1  exp(-1i*regul*(spacing(1)/sum(spacing))) exp(-2*1i*regul*(spacing(1)/sum(spacing)));
        exp(2*1i*regul*((spacing(1)+spacing(2))/sum(spacing))) exp(1i*regul*((spacing(1)+spacing(2))/sum(spacing))) 1  exp(-1i*regul*((spacing(1)+spacing(2))/sum(spacing))) exp(-2i*regul*((spacing(1)+spacing(2))/sum(spacing)));
        exp(2*1i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing))) exp(1i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing))) 1  exp(-1i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing))) exp(-2i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing)));
        exp(2*1i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing))) exp(1i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing))) 1  exp(-1i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing))) exp(-2i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing)));];%488
    phase_matrix = inv(phase_matrix);
    ztoxy = zeros(size(phase_im,4),size(phase_im,1)*size(phase_im,2)*size(phase_im,3));
    sep_im = zeros(size(phase_im));
    for testi=1:size(phase_im,4)       
        F = phase_im(:,:,:,testi);
        ztoxy(testi,:) = F(:);
    end
    [sizex, sizey, sizez, nangle] = size(phase_im);
    clear F phase_im
    ztoxy = phase_matrix * ztoxy; 
    for testi=1:nangle          
        sep_im(:,:,:,testi) = reshape(ztoxy(testi,:),[sizex, sizey, sizez]); 
        sep_im(:,:,:,testi) = fftshift(fftn(sep_im(:,:,:,testi)));
    end
end