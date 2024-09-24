% Cross correlation Method:GPU accelerate Pattern wave vector,p0 OTF and p1 OTF in 2D
function [zuobiaotmpx,zuobiaotmpy,tmpangle6,tmpc6,EvaluateParam,Kurtosis]=SIM_3D_p_CC_saveHistg(spzhongxin, spyiweiP, spzhongxinAC, spyiweiPAC, spyiweiPhalfAC,...
            H_zhongxin, H_yiweiP, H_yiweiPhalf,fc, pg, sizex, sizey, fanwei, orti, jindu)
hangs=2*orti-1;
n_512 = max([sizex,sizey]);
n = 512;
pg = pg*(n_512/n);
fc = fc*(n_512/n);
Pindex = 1;
Phalf_index = 2;
[otfx, otfy] = size(H_zhongxin);
zuobiaotmpx(1:9,1:2) = n_512;              
zuobiaotmpy(1:9,1:2) = n_512;               
[k_x_512, k_y_512] = meshgrid(-(n_512)/2:(n_512)/2-1, -(n_512)/2:(n_512)/2-1);
k_r_512 = sqrt(k_x_512.^2+k_y_512.^2);    
indi_512 =  k_r_512 > fc ;              
K_h = [otfx, otfy];
N_h = 2*K_h;
L_h = ceil((N_h-K_h) / 2);
v_h = colonvec(L_h+1, L_h+K_h);
hw = zeros(N_h);
H_zhongxin(indi_512) = 0;                  
H_yiweiP(indi_512) = 0;
H_yiweiPhalf(indi_512) = 0;
H_zhongxin = abs(H_zhongxin);
H_yiweiP = abs(H_yiweiP);
H_yiweiPhalf = abs(H_yiweiPhalf);
% m=0 OTF
H1_zhongxin = H_zhongxin;                  
H1_zhongxin(H1_zhongxin~=0) = 1;
hw(v_h{:}) = H_zhongxin;                   
H = hw;
hw(v_h{:}) = H1_zhongxin;
H1 = hw;
H = gpuArray(single(H));
H1 = gpuArray(single(H1));
% m=2 OTF
H1_yiweiP = H_yiweiP;                        
H1_yiweiP(H1_yiweiP~=0) = 1;
hw(v_h{:}) = H_yiweiP;                    
H_yiweiP = hw;
hw(v_h{:}) = H1_yiweiP;
H1_yiweiP = hw;
H_yiweiP = gpuArray(single(H_yiweiP));
H1_yiweiP = gpuArray(single(H1_yiweiP));
% m=1 OTF
H1_yiweiPhalf = H_yiweiPhalf;             
H1_yiweiPhalf(H1_yiweiPhalf~=0) = 1;
hw(v_h{:}) = H_yiweiPhalf;                 
H_yiweiPhalf = hw;
hw(v_h{:}) = H1_yiweiPhalf;
H1_yiweiPhalf = hw;
H_yiweiPhalf = gpuArray(single(H_yiweiPhalf));
H1_yiweiPhalf = gpuArray(single(H1_yiweiPhalf));

%% calculation
buchangx = 1;
buchangy = 1;
Hzuan = ones(sizex, sizey);
cishu = dft(Hzuan,Hzuan);               
sp_tmp = conj(spyiweiP);
sp_tmp = sp_tmp(end:-1:1,end:-1:1);      
jieguo = dft(spzhongxin,sp_tmp);         
lihe = abs(jieguo./(cishu+eps));
[k_x_1024, k_y_1024] = meshgrid(-(2*n_512)/2+1:(2*n_512)/2-1, -(2*n_512)/2+1:(2*n_512)/2-1);
k_r_1024 = sqrt(k_x_1024.^2+k_y_1024.^2);
indi = ((k_r_1024 < (floor(pg/Pindex)-fanwei))|(k_r_1024 > (floor(pg/Pindex)+fanwei)));  
lihe(indi) = 0;
[xx,yy] = ind2sub(size(lihe),find(lihe==max(lihe(:)))); 
EvaluateParam.CCMax = max(lihe(:));
% Gradient descent
maxx = xx;
maxy = yy;
ky = 2*pi*(maxx-n_512)/(2*n_512);
kx = 2*pi*(maxy-n_512)/(2*n_512);
[xx2,yy2] = meshgrid(-n_512:n_512-1, -n_512:n_512-1);  
xx2 = gpuArray(single(xx2));
yy2 = gpuArray(single(yy2));
% if orti == 3 || orti == 2
%     imwritestacka(lihe, 'nihe.tif');
% elseif orti == 1
%     imwritestack(lihe, 'nihe.tif');
% end
hw(v_h{:}) = spyiweiP;    
spyiweiP = hw;
spyiweiP = gpuArray(single(spyiweiP));
hw(v_h{:}) = spzhongxin;    
spzhongxin = hw;
spzhongxin = gpuArray(single(spzhongxin));
Irtest = exp(1i*(kx*xx2+ky*yy2));
replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiP))).*Irtest))); 
replcHtest = abs(replcHtest);
replcHtest(abs(replcHtest)>0.9) = 1; 
replcHtest(abs(replcHtest)~=1) = 0;
replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiP))).*Irtest))));  
replch = replch.*replcHtest; 
hetest = (conj(fftshift(fft2(ifft2(ifftshift(spyiweiP)).*Irtest)).*replcHtest.*H)).*(spzhongxin.*replch);
hetest = abs(sum(hetest(:)));
cishuhtest = H1.*replcHtest; 
cishuhtest = sum(cishuhtest(:));
hetest = hetest./cishuhtest;
he = hetest;
% Determining the direction of the gradient
maxx_tmp1 = maxx-10^-5;
maxx_tmp2 = maxx+10^-5;
maxy_tmp1 = maxy-10^-5;
maxy_tmp2 = maxy+10^-5;
ky_tmp1 = 2*pi*(maxx_tmp1-n_512)/(2*n_512);
kx_tmp1 = 2*pi*(maxy_tmp1-n_512)/(2*n_512);
ky_tmp2 = 2*pi*(maxx_tmp2-n_512)/(2*n_512);
kx_tmp2 = 2*pi*(maxy_tmp2-n_512)/(2*n_512);
for ii=1:1:4  % Finding the optimum in different directions
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
    replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiP))).*Irtest)));
    replcHtest = abs(replcHtest);
    replcHtest(abs(replcHtest)>0.9) = 1;
    replcHtest(abs(replcHtest)~=1) = 0;
    replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiP))).*Irtest))));
    replch = replch.*replcHtest;
    hetest = (conj(fftshift(fft2(ifft2(ifftshift(spyiweiP)).*Irtest)).*replcHtest.*H)).*(spzhongxin.*replch);
    hetest = abs(sum(hetest(:)));
    cishuhtest = H1.*replcHtest;
    cishuhtest = sum(cishuhtest(:));
    hetest = hetest./cishuhtest;
    test(ii) = hetest;
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
while ((buchangx>jindu)||(buchangy>jindu))
    %%%%%%%%%%%%%%%%%%%%%%%maxxorientation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxx_tmp1 = maxx-10^-5;
    maxx_tmp2 = maxx+10^-5;
    ky_tmp1 = 2*pi*(maxx_tmp1-n_512)/(2*n_512);
    ky_tmp2 = 2*pi*(maxx_tmp2-n_512)/(2*n_512);
    for ii=3:1:4   % x orientation £º-1 +1
        switch ii
            case 3
                kxtest = 2*pi*(maxy-n_512)/(2*n_512);
                kytest = ky_tmp1;
            case 4
                kxtest = 2*pi*(maxy-n_512)/(2*n_512);
                kytest = ky_tmp2;
        end
        Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
        replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiP))).*Irtest))); 
        replcHtest = abs(replcHtest);
        replcHtest(abs(replcHtest)>0.9) = 1;
        replcHtest(abs(replcHtest)~=1) = 0;
        replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiP))).*Irtest)))); 
        replch = replch.*replcHtest;  
        hetest = (conj(fftshift(fft2(ifft2(ifftshift(spyiweiP)).*Irtest)).*replcHtest.*H)).*(spzhongxin.*replch); 
        hetest = abs(sum(hetest(:)));
        cishuhtest = H1.*replcHtest;
        cishuhtest = sum(cishuhtest(:));
        hetest = hetest./cishuhtest;
        test(ii) = hetest;
    end
    if((test(3)>test(4)))
        flag_maxx=-1;
    elseif((test(3)<test(4)))
        flag_maxx=+1;
    else
        flag_maxx=-1* flag_maxx;
    end
    while(buchangx>(jindu))
        maxx_tmp = maxx+flag_maxx*buchangx;
        kytest = 2*pi*(maxx_tmp-n_512)/(2*n_512);
        kxtest = 2*pi*(maxy-n_512)/(2*n_512);
        Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
        replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiP))).*Irtest)));
        replcHtest = abs(replcHtest);
        replcHtest(abs(replcHtest)>0.9)=1;
        replcHtest(abs(replcHtest)~=1)=0;
        replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiP))).*Irtest))));
        replch = replch.*replcHtest;
        hetest = (conj(fftshift(fft2(ifft2(ifftshift(spyiweiP)).*Irtest)).*replcHtest.*H)).*(spzhongxin.*replch);
        hetest = abs(sum(hetest(:)));
        cishuhtest = H1.*replcHtest;
        cishuhtest = sum(cishuhtest(:));
        hetest = hetest./cishuhtest;
        he_tmp = hetest;
        if(he_tmp<=he)
            buchangx=0.5*buchangx;
        elseif(he_tmp>he)
            he=he_tmp;      
            maxx=maxx_tmp;  
            break;
        end
    end    % end of while x
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%maxyorientation%%%%%%%%%%%%%%%%%%%%
    maxy_tmp1=maxy-10^-5;
    maxy_tmp2=maxy+10^-5;
    kx_tmp1 = 2*pi*(maxy_tmp1-n_512)/(2*n_512);
    kx_tmp2 = 2*pi*(maxy_tmp2-n_512)/(2*n_512);
    for ii=1:1:2   % y orientation£º-1 +1
        switch ii
            case 1
                kxtest = kx_tmp1;
                kytest = 2*pi*(maxx-n_512)/(2*n_512);
            case 2
                kxtest = kx_tmp2;
                kytest = 2*pi*(maxx-n_512)/(2*n_512);
        end
        Irtest=exp(1i*(kxtest*xx2+kytest*yy2));
        replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiP))).*Irtest)));
        replcHtest = abs(replcHtest);
        replcHtest(abs(replcHtest)>0.9)=1;
        replcHtest(abs(replcHtest)~=1)=0;
        replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiP))).*Irtest))));
        replch= replch.*replcHtest;
        hetest= (conj(fftshift(fft2(ifft2(ifftshift(spyiweiP)).*Irtest)).*replcHtest.*H)).*(spzhongxin.*replch);
        hetest=abs(sum(hetest(:)));
        cishuhtest=H1.*replcHtest;
        cishuhtest=sum(cishuhtest(:));
        hetest=hetest./cishuhtest;
        test(ii)=hetest;
    end
    if((test(1)>test(2)))
        flag_maxy=-1;
    elseif((test(1)<test(2)))
        flag_maxy=+1;
    else
        flag_maxy=-1*flag_maxy;
    end
    while(buchangy>(jindu))
        maxy_tmp = maxy+flag_maxy*buchangy;
        kytest = 2*pi*(maxx-n_512)/(2*n_512);
        kxtest = 2*pi*(maxy_tmp-n_512)/(2*n_512);
        Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
        replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiP))).*Irtest)));
        replcHtest = abs(replcHtest);
        replcHtest(abs(replcHtest)>0.9) = 1;
        replcHtest(abs(replcHtest)~=1)=0;
        replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiP))).*Irtest))));
        replch = replch.*replcHtest;
        hetest = (conj(fftshift(fft2(ifft2(ifftshift(spyiweiP)).*Irtest)).*replcHtest.*H)).*(spzhongxin.*replch);
        hetest = abs(sum(hetest(:)));
        cishuhtest = H1.*replcHtest;
        cishuhtest = sum(cishuhtest(:));
        hetest = hetest./cishuhtest;
        he_tmp = hetest;
        if(he_tmp<=he)
            buchangy = 0.5*buchangy;
        elseif(he_tmp>he)
            he = he_tmp;    
            maxy = maxy_tmp;
            break;
        end
    end       % end of while y
end
hw(v_h{:}) = spyiweiPAC;       
spyiweiPAC = hw;
hw(v_h{:}) = spzhongxinAC;   
spzhongxinAC = hw;
hw(v_h{:}) = spyiweiPhalfAC;   
spyiweiPhalfAC = hw;
clearvars -except spyiweiPAC Irtest replcHtest H maxx maxy Pindex hangs zuobiaotmpx zuobiaotmpy n_512 tmpangle6 tmpc6 ...
    spzhongxinAC replch youhuatest orti xx2 yy2 Phalf_index spyiweiPhalfAC H1_yiweiPhalf H_yiweiPhalf isRealTime EvaluateParam

% m=2 parameters
youhuatest = gather(fftshift(fft2(ifft2(ifftshift(spyiweiPAC)).*Irtest)).*replcHtest.*H);
zuobiaotmpx(hangs+1+floor((hangs-1)/2), Pindex) = maxx;           
zuobiaotmpy(hangs+1+floor((hangs-1)/2), Pindex) = maxy;
zuobiaotmpx(hangs+2+floor((hangs-1)/2), Pindex) = 2*n_512-maxx;   
zuobiaotmpy(hangs+2+floor((hangs-1)/2), Pindex) = 2*n_512-maxy;
[tmpangle6(Pindex), tmpc6(Pindex), Kurtosis(:,1)] = SIM_3D_angle_EMD(gather(spzhongxinAC),gather(replch),youhuatest,gather(H),orti,Pindex);

% m=1 parameters
maxx = n_512 + (maxx-n_512)/2;
maxy = n_512 + (maxy-n_512)/2;
kytest = 2*pi*(maxx-n_512)/(2*n_512);
kxtest = 2*pi*(maxy-n_512)/(2*n_512);
Irtest = exp(1i*(kxtest*xx2+kytest*yy2));
replcHtest = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_yiweiPhalf))).*Irtest)));
replcHtest = abs(replcHtest);
replcHtest(abs(replcHtest)>0.9) = 1;
replcHtest(abs(replcHtest)~=1) = 0;
replch = abs(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_yiweiPhalf))).*Irtest))));
replch = replch.*replcHtest;
youhuatest = gather(fftshift(fft2(ifft2(ifftshift(spyiweiPhalfAC)).*Irtest)).*replcHtest.*H);
zuobiaotmpx(hangs+1+floor((hangs-1)/2),Phalf_index) = maxx;          
zuobiaotmpy(hangs+1+floor((hangs-1)/2),Phalf_index) = maxy;
zuobiaotmpx(hangs+2+floor((hangs-1)/2),Phalf_index) = 2*n_512-maxx;   
zuobiaotmpy(hangs+2+floor((hangs-1)/2),Phalf_index) = 2*n_512-maxy;

[tmpangle6(Phalf_index),tmpc6(Phalf_index),Kurtosis(:,2)] = SIM_3D_angle_EMD(gather(spzhongxinAC),gather(replch),youhuatest,gather(H),orti,Phalf_index);
