% Cross correlation Method:GPU accelerate Pattern wave vector,p0 OTF and p1 OTF in 2D
function [tmpangle6,tmpc6,Kurtosisa]=SIM_3D_p_CCa_saveHistg(angle,H_zhongxin,fc, sizex, sizey, orti, spzhongxina,spyiweiPhalfa,H_zhongxina,H_yiweiPhalfa,zuobiaox,zuobiaoy)
hangs=2*orti-1;
n_512 = max([sizex,sizey]);
Pindex = 1;
Phalf_index = 2;
[otfx, otfy] = size(H_zhongxin);
% zuobiaotmpx(1:9,1:2) = n_512;             
% zuobiaotmpy(1:9,1:2) = n_512;               
[k_x_512, k_y_512] = meshgrid(-(n_512)/2:(n_512)/2-1, -(n_512)/2:(n_512)/2-1);
k_r_512 = sqrt(k_x_512.^2+k_y_512.^2);     
indi_512 =  k_r_512 > fc ;             
K_h = [otfx, otfy];
N_h = 2*K_h;
L_h = ceil((N_h-K_h) / 2);
v_h = colonvec(L_h+1, L_h+K_h);
hw = zeros(N_h);

% pa 
H_yiweiPhalfa(indi_512)=0;
H_zhongxina(indi_512)=0;
H_yiweiPhalfa=abs(H_yiweiPhalfa);
H_zhongxina=abs(H_zhongxina);

% pa m=1 OTF
H1_yiweiPhalfa = H_yiweiPhalfa;              
H1_yiweiPhalfa(H1_yiweiPhalfa~=0) = 1;
hw(v_h{:}) = H_yiweiPhalfa;             
H_yiweiPhalfa = hw;
hw(v_h{:}) = H1_yiweiPhalfa;
H1_yiweiPhalfa = hw;
H_yiweiPhalfa = gpuArray(single(H_yiweiPhalfa));
H1_yiweiPhalfa = gpuArray(single(H1_yiweiPhalfa));

% pa m=0 OTF
H1_zhongxina = H_zhongxina;                  
H1_zhongxina(H1_zhongxina~=0) = 1;
hw(v_h{:}) = H_zhongxina;                  
Ha = hw;
hw(v_h{:}) = H1_zhongxina;
H1a = hw;
Ha = gpuArray(single(Ha));
H1a = gpuArray(single(H1a));

%% calculation
hw(v_h{:}) = spzhongxina;    
spzhongxina = hw;
hw(v_h{:}) = spyiweiPhalfa;   
spyiweiPhalfa = hw;

%% m=1 parameters
maxx = zuobiaox(hangs+1+floor((hangs-1)/2), Pindex);
maxy = zuobiaoy(hangs+1+floor((hangs-1)/2), Pindex);
[xx2,yy2] = meshgrid(-n_512:n_512-1, -n_512:n_512-1);
xx2 = gpuArray(single(xx2));
yy2 = gpuArray(single(yy2));
maxx = n_512 + (maxx-n_512)/2;
maxy = n_512 + (maxy-n_512)/2;
kytest = 2*pi*(maxx-n_512)/(2*n_512);
kxtest = 2*pi*(maxy-n_512)/(2*n_512);
Irtest = exp(1i*(kxtest*xx2+kytest*yy2));

% pa estimation
phalfa=exp(-1i*angle);
replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiPhalfa))).*Irtest)));
replcHtest = abs(replcHtest);
replcHtest(abs(replcHtest)>0.9) = 1;
replcHtest(abs(replcHtest)~=1) = 0;
replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiPhalfa))).*Irtest))));
replch = replch.*replcHtest;
youhuatest = gather(fftshift(fft2(ifft2(ifftshift(spyiweiPhalfa)).*Irtest)).*replcHtest.*Ha.*phalfa);

[tmpangle6(Phalf_index),tmpc6(Phalf_index),Kurtosisa(:,1)] = SIM_3D_angle_EMD(gather(spzhongxina),gather(replch),youhuatest,gather(Ha),orti,Phalf_index);


