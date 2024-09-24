function SIM_4PS_Wiener_deconv(H_all,pathname,filename,sizex,sizey,numz,nort,nangle,...
    bgflag,bgname,zuobiaox,zuobiaoy,angle6,c6,regul,spjg,weilac,notch_swith,notch_para_1,notch_para_2,fc_bhs,...
    ~,pixelsizeZ,wavelengh,Zpadding,isEqualization,I2M_para)
% OTF
H_all = single(H_all);
H = H_all(:,:,1:numz);            % m=0 OTF
hp = H_all(:,:,2*numz+1:numz*3);  % m=1 OTF
h2p = H_all(:,:,1*numz+1:numz*2);    % m=2 OTF
clear H_all
n = size(H,1);
% resize OTF
H = imresize3(H, [n, n, numz+2*Zpadding], 'linear');
hp = imresize3(hp, [n, n, numz+2*Zpadding], 'linear');
h2p = imresize3(h2p, [n, n, numz+2*Zpadding], 'linear');

all_name = [pathname,filename];

% background noise
if bgflag == 1
    bg = single(imreadstack(bgname));
    [bgsizex, bgsizey] = size(bg);
    if bgsizex == sizex && bgsizey == sizey
    else
        bg = ones(sizex, sizey, 'single').*100;
    end
elseif bgflag == 0
    bg = zeros(sizex,sizey, 'single');
end

% Scaling coordinate
zuobiaox = zuobiaox*(n/zuobiaox(1,1));
zuobiaoy = zuobiaoy*(n/zuobiaoy(1,1));

% read parameters
xishu = [(c6(1,1));(c6(1,2));1;(c6(1,2));(c6(1,1));
    (c6(2,1));(c6(2,2));1;(c6(2,2));(c6(2,1));
    (c6(3,1));(c6(3,2));1;(c6(3,2));(c6(3,1));];
fc_bhs = ceil(fc_bhs*(n/512));  % maximum spatial frequency
plong = sum(((zuobiaox-zuobiaox(1,1)).^2+(zuobiaoy-zuobiaoy(1,1)).^2).^0.5, 1)./6;
deph1 = angle6(1,1);    % m=2 initial phase
deph2 = angle6(2,1);
deph3 = angle6(3,1);
deph4 = angle6(1,2);    % m=1 initial phase
deph5 = angle6(2,2);
deph6 = angle6(3,2);
reh = zeros((2*n), (2*n), (numz+2*Zpadding), 'single');   
Irtest = zeros((2*n), (2*n), nort*nangle, 'single');      
[xx2, yy2] = meshgrid(-n:n-1, -n:n-1);                    
[k_x, k_y] = meshgrid(-(n)/2:(n)/2-1, -(n)/2:(n)/2-1);
k_r = sqrt(k_x.^2+k_y.^2);
clear k_x k_y

% Transform the order of coordinates
[tmp_zuobiaox,tmp_zuobiaoy] = zub_change(zuobiaox,zuobiaoy);  

% OTF and OTF mask
% Mask_4PS=Make_4PS_mask(); % mask with 8 um depth and 40 nm step size
Mask_4PS=imreadstack('Data/OTF/Mask_4PS.tif');
Mask_4PS = GenerateOTF_Fre_revised(Mask_4PS, pixelsizeZ, pixelsizeZ, numz, n);
otfmask0 = Mask_4PS(:,:,1:size(Mask_4PS,3)/3);
otfmask2 = Mask_4PS(:,:,(size(Mask_4PS,3)/3+1):2*size(Mask_4PS,3)/3);
otfmask1 = Mask_4PS(:,:,(2*size(Mask_4PS,3)/3+1):size(Mask_4PS,3));
clear Mask_4PS

otfmask0 = imresize3(single(otfmask0), [n, n, numz+2*Zpadding], 'linear');
otfmask0(otfmask0<0.1) = 0;
otfmask0(otfmask0~=0) = 1;

otfmask1 = imresize3(single(otfmask1), [n, n, numz+2*Zpadding], 'linear');
otfmask1(otfmask1<0.1) = 0;
otfmask1(otfmask1~=0) = 1;

otfmask2 = imresize3(single(otfmask2), [n, n, numz+2*Zpadding], 'linear');
otfmask2(otfmask2<0.1) = 0;
otfmask2(otfmask2~=0) = 1;

Hk = H.*otfmask0;          % m=0 beyond the cut-off frequency is 0
h2pPad = h2p.*otfmask2;    % m=2
hpPad = hp.*otfmask1;      % m=1

% zero padding
K_h = size(Hk);
if numel(K_h) == 3
    N_h = 2*K_h-[0,0,K_h(1,3)];
else   
    N_h = 2*K_h-[0,0];
end
L_h = ceil((N_h-K_h) / 2);
v_h = colonvec(L_h+1, L_h+K_h);
hw = zeros(N_h, 'single');
hw(v_h{:}) = Hk;
Hk = single(hw);
hw(v_h{:}) = h2pPad;
h2pPad = single(hw);
hw(v_h{:}) = hpPad;
hpPad = single(hw);

clear hw N_h K_h L_h v_h

% notch filter
if notch_swith==1
    notch_filter = 1-exp(-notch_para_1*(abs(k_r).^notch_para_2));
    notch_filter = single(repmat(notch_filter,[1,1,nangle]));
elseif notch_swith==0
    notch_filter = ones(size(k_r), 'single');
    notch_filter = repmat(notch_filter,[1,1,nangle]);
end
% imwritestack(notch_filter, './konyu_notch_filter.tif');
notch_filter = single(notch_filter);

% OTF mask binarized area
for ii=1:(nort*nangle)
    kytest = 2*pi*(tmp_zuobiaox(ii)-n)/(2*n);
    kxtest = 2*pi*(tmp_zuobiaoy(ii)-n)/(2*n);
    Irtest(:,:,ii) = exp(1i*(kxtest*xx2+kytest*yy2)); 
end

% apodization function
[bhs] = bhs_4PS_theory(n, numz, plong, fc_bhs,pixelsizeZ,wavelengh,I2M_para); 
% imwritestack(bhs,'bhs.tif')
bhs = single(imresize3(bhs, [2*n, 2*n, numz+2*Zpadding], 'linear'));

% remove background and OTF shift
disp('Calculating Wiener denominator, please wait...');
hs = zeros((2*n), (2*n), (numz+2*Zpadding),  'single');
for ii = 1:(nort*nangle)
    if mod(ii,nangle)==ceil(nangle/2)-1 || mod(ii,nangle)==ceil(nangle/2)+1 
        reh = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(xishu(ii).*hpPad))).*Irtest(:,:,ii)) ));       % OTF  m=1
        OTFmask=otfmask1;
    elseif mod(ii,nangle)==ceil(nangle/2)
        reh = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(xishu(ii).*Hk))).*Irtest(:,:,ii)) ));          % OTF  m=0
        OTFmask=otfmask0;
    else
        reh = fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(xishu(ii).*h2pPad))).*Irtest(:,:,ii)) ));      % OTF  m=2
        OTFmask=otfmask2;
    end

    OTFmask = imresize3(single(OTFmask), [2*n, 2*n, numz+2*Zpadding], 'linear');
    OTFmask=fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift(OTFmask))).*Irtest(:,:,ii)) ));
    OTFmask(OTFmask<0.1) = 0;
    OTFmask(OTFmask~=0) = 1;

    reh = abs(reh);
    reh = reh.*OTFmask; 
    hs = hs + reh.^2;
end
clear reh OTFmask qx qy qz qpar otfmask1 otfmask0 otfmask2
hs_Shift = zeros(n, n, (numz+2*Zpadding), (nort*nangle), 'single');
bhs_Shift = zeros(n, n, (numz+2*Zpadding), (nort*nangle), 'single');
for ii= 1:(nort*nangle)
    hs_Shift(:,:,:,ii) = Crop(fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift( hs ))).*conj(Irtest(:,:,ii))))),H);
    bhs_Shift(:,:,:,ii) = Crop(fftshift(fftn(ifftshift(fftshift(ifftn(ifftshift( bhs ))).*conj(Irtest(:,:,ii))))),H);
end
WienerDenom_Move = abs(bhs_Shift)./(abs(hs_Shift) + 0.005*(weilac)^2);
clear hs hs_Shift bhs_Shift mask sigma x y

% Separation phase matrix
phase_matrix= [1 1 1 1 1;
    exp(2*1i*regul*(spjg(1)/sum(spjg)))  exp(1i*regul*(spjg(1)/sum(spjg))) 1  exp(-1i*regul*(spjg(1)/sum(spjg))) exp(-2*1i*regul*(spjg(1)/sum(spjg)));
    exp(2*1i*regul*((spjg(1)+spjg(2))/sum(spjg))) exp(1i*regul*((spjg(1)+spjg(2))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2))/sum(spjg))) exp(-2i*regul*((spjg(1)+spjg(2))/sum(spjg)));
    exp(2*1i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg))) exp(1i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg))) exp(-2i*regul*((spjg(1)+spjg(2)+spjg(3))/sum(spjg)));
    exp(2*1i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg))) exp(1i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg))) exp(-2i*regul*((spjg(1)+spjg(2)+spjg(3)+spjg(4))/sum(spjg)));];
phase_matrix1=phase_matrix.*repmat([exp(1i*deph1),exp(1i*deph4),1,exp(-1i*deph4),exp(-1i*deph1)],[size(phase_matrix,1),1,1]);
phase_matrix2=phase_matrix.*repmat([exp(1i*deph2),exp(1i*deph5),1,exp(-1i*deph5),exp(-1i*deph2)],[size(phase_matrix,1),1,1]);
phase_matrix3=phase_matrix.*repmat([exp(1i*deph3),exp(1i*deph6),1,exp(-1i*deph6),exp(-1i*deph3)],[size(phase_matrix,1),1,1]);
phase_matrix1 = inv(phase_matrix1);
phase_matrix2 = inv(phase_matrix2);
phase_matrix3 = inv(phase_matrix3);
phase_matrix(:,:,1) = phase_matrix1;
phase_matrix(:,:,2) = phase_matrix2;
phase_matrix(:,:,3) = phase_matrix3;
padsize = 0;
x = 1:(sizey+2*padsize);
y = (1:(sizex+2*padsize))';
sigma = 0.25;
mask = single(repmat(sigmoid(sigma*(x-padsize)) - sigmoid(sigma*(x-sizey-padsize-1)), sizex+2*padsize, 1) .* repmat(sigmoid(sigma*(y-padsize)) - sigmoid(sigma*(y-sizex-padsize-1)), 1, sizey+2*padsize));
K_h2 = [n, n, (numz+2*Zpadding)];
N_h2 = [2*n, 2*n, (numz+2*Zpadding)];
L_h2 = ceil((N_h2-K_h2) / 2);
v_h2 = colonvec(L_h2+1, L_h2+K_h2);
K_h1 = [sizex,sizey,(numz+2*Zpadding)];
N_h1 = [n,n,(numz+2*Zpadding)];
L_h1 = ceil((N_h1-K_h1) / 2);
v_h1 = colonvec(L_h1+1, L_h1+K_h1);
hw1 = zeros(N_h1, 'single');

clearvars -except nort nangle numz all_name pathname filename sizex sizey bg notch_filter mask hs ...
    Irtest bhs n v_h phase_matrix v_h2 xishu N_h2 H hp WienerDenom_Move jiequ Zpadding weilac hw1 v_h1 h2pPad h2p isEqualization

num = 0;
num_images = nort*nangle*numz; 
wd = Window(numz+2*Zpadding, sizex, sizey, 0.25, 3);
hw1(v_h1{:}) = wd;
wd = hw1;
clear hw1

spsim = zeros(n, n, (numz+2*Zpadding), nort*nangle, 'single');
fimage = zeros(2*n, 2*n, (numz+2*Zpadding), 'single');
num = num + 1;
if ( 1+(num-1)*(nort*nangle*numz)+(nort*nangle*numz)-1 <= num_images )
    SIM_raw = single(myimreadstack_TIRF(all_name, 1+(num-1)*(nort*nangle*numz), (nort*nangle*numz), sizex, sizey));
else
    error('Length overflows');
end
SIM_raw = SIM_raw - bg;
SIM_raw(SIM_raw<0) = 0;

% Balanced image
if isEqualization == 1
    disp('Performing equalization in stacks, please wait...');
    restack=zeros(sizex,sizey,numz,nangle*nort);
    weights=zeros(1,nangle*nort);
    for ii=1:nort
        for jj=1:nangle
            restack(:,:,1:numz,(ii-1)*nangle+jj) = SIM_raw(:,:,(ii-1)*nangle+jj:nangle*nort:numz*nangle*nort);
            weights(:,(ii-1)*nangle+jj) = mean(mean(mean(restack(:,:,:,(ii-1)*nangle+jj))));
        end
    end
    clear restack
    weights = mean(weights)./weights;
    disp(['Balanced weight: ' num2str(weights)]);
end
SIM_raw = SIM_raw.*(mask.^3);

% Separation in spatial domain
disp('Seperating raw data, please wait...');
for orti=1:nort
    phase_im = zeros(size(SIM_raw,1),size(SIM_raw,2),(numz+2*Zpadding),nangle,'single');
    for anglei=1:nangle 
        if isEqualization == 0
            phase_im(:,:,Zpadding+1:end-Zpadding,anglei) = SIM_raw(:,:,anglei+(orti-1)*nangle:nangle*nort:numz*nangle*nort);
        else
            phase_im(:,:,Zpadding+1:end-Zpadding,anglei) = SIM_raw(:,:,anglei+(orti-1)*nangle:nangle*nort:numz*nangle*nort).*weights(:,anglei+(orti-1)*nangle);
        end

        if Zpadding~=0
            phase_im(:,:,1:Zpadding,anglei) = repmat(SIM_raw(:,:,anglei+(orti-1)*nangle),[1,1,Zpadding]);
            phase_im(:,:,end-Zpadding+1:end,anglei) = repmat(SIM_raw(:,:,numz*nangle*nort-(nangle*nort-(anglei+(orti-1)*nangle))),[1,1,Zpadding]);
        end
        phase_im(:,:,:,anglei) = phase_im(:,:,:,anglei).*wd;
    end
    ztoxy = zeros(size(phase_im,4),size(phase_im,1)*size(phase_im,2)*size(phase_im,3),'single');
    for testi=1:size(phase_im,4)
        F = (phase_im(:,:,:,testi));
        ztoxy(testi,:) = F(:);
    end
    clear F
    ztoxy = phase_matrix(:,:,orti) * ztoxy;
    for testi=1:size(phase_im,4)  
        spsim(:,:,:,(orti-1)*nangle + testi) = reshape(ztoxy(testi,:),[size(phase_im,1),size(phase_im,2),size(phase_im,3)]); %将每一帧数据存为一行，n帧三维数据即转化为n行二维进行处理
        spsim(:,:,:,(orti-1)*nangle + testi) = fftshift(fftn(spsim(:,:,:,(orti-1)*nangle + testi)));
        spsim(:,:,:,(orti-1)*nangle + testi) = spsim(:,:,:,(orti-1)*nangle + testi).* notch_filter(:,:,testi);
    end
end
clear ztoxy phase_im SIM_raw sep_im

% recombine filtered bands
disp('Performing Wiener filtering, please wait...');
hw2 = zeros(N_h2,'single');
for t = 1:nort*nangle
    if mod(t,nangle) == ceil(nangle/2)-1 || mod(t,nangle)==ceil(nangle/2)+1
        spsim(:,:,:,t) = xishu(t).* spsim(:,:,:,t).* conj(hp).* WienerDenom_Move(:,:,:,t);
    elseif mod(t,nangle) == ceil(nangle/2)
        spsim(:,:,:,t) = xishu(t).* spsim(:,:,:,t).* conj(H).* WienerDenom_Move(:,:,:,t);
    else
        spsim(:,:,:,t) = xishu(t).* spsim(:,:,:,t).* conj(h2p).* WienerDenom_Move(:,:,:,t);
    end
    hw2(v_h2{:}) = spsim(:,:,:,t);
    fimage = fimage + ifftn(ifftshift(hw2)).*Irtest(:, :, t);
end
clear spsim hw2
fimage = real(fimage);
fimage = fimage((end/2)+1-(sizex):(end/2)+(sizex),(end/2)+1-(sizey):(end/2)+(sizey),:);
fimage(fimage<0)=0;

% save data
if num == 1
    imwritestack(abs( fftshift( fftn(fimage(:,:,Zpadding+1:end-Zpadding)) ) ), [pathname 'SIM Result/SpectrumWiener' filename(1:end-4) '.tif']);
else
    imwritestacka(abs( fftshift( fftn(fimage(:,:,Zpadding+1:end-Zpadding)) )), [pathname 'SIM Result/SpectrumWiener' filename(1:end-4) '.tif']);
end
if num == 1
    imwritestack(fimage(:,:,Zpadding+1:end-Zpadding), [pathname 'SIM Result/Wiener' filename(1:end-4) '.tif']);
else
    imwritestacka(fimage(:,:,Zpadding+1:end-Zpadding), [pathname 'SIM Result/Wiener' filename(1:end-4) '.tif']);
end
clear fimage
disp('The 4Pi-SIM image reconstruction is completed.');
end