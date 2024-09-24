function  [sep_im]=SIM_3D_seprt_PAZ(SIM_raw,nangle,nort,spjg,regul,numz,orti)
n = max([size(SIM_raw,1),size(SIM_raw,2)]);
K_h = [size(SIM_raw,1),size(SIM_raw,2)];
N_h = [n, n];
L_h = ceil((N_h-K_h) / 2);
v_h = colonvec(L_h+1, L_h+K_h);
hw = zeros(N_h);
fd512 = zeros(n, n, size(SIM_raw,3));
for ii=1:1:size(SIM_raw,3)  
    hw(v_h{:}) = SIM_raw(:,:,ii);
    fd512(:,:,ii) = hw;
end
SIM_raw = fd512;
clear fd512 hw;
phase_im = zeros(size(SIM_raw,1), size(SIM_raw,2), numz, nangle);
for anglei=1:nangle  
    phase_im(:,:,:,anglei) = SIM_raw(:,:,anglei+(orti-1)*nangle:nangle*nort:numz*nangle*nort);
end
clear SIM_raw
phase_matrix= [1 1 1 1 1; % separation phase matrix
    exp(2*1i*regul*(spjg(1)/sum(spjg)))  exp(1i*regul*(spjg(1)/sum(spjg))) 1  exp(-1i*regul*(spjg(1)/sum(spjg))) exp(-2*1i*regul*(spjg(1)/sum(spjg)));
    exp(2*1i*regul*((spjg(1)+spjg(2))/sum(spjg))) exp(1i*regul*((spjg(1)+spjg(2))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2))/sum(spjg))) exp(-2i*regul*((spjg(1)+spjg(2))/sum(spjg)));
    exp(2*1i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg))) exp(1i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg))) exp(-2i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg)));
    exp(2*1i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg))) exp(1i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg))) exp(-2i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg)));];%488
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