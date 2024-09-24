% mask 4Pi-SIM OTF mask 8um 40nm step
function Mask_4PS=Make_4PS_mask()
pixelsize=63.7255;
pixelsizeZ=40;
lambdaex=488;
lambda=525;
Pzrefmed=1.333;
n=512;
numz=200;
Rawpixelsize = [pixelsize, pixelsize, pixelsizeZ];
OTFmaskrefmed=1.27;
ReceiveNA=1.1;
xmedian=floor(n/2)+1;
zmedian=floor(numz/2)+1;

% I2M OTF mask
DqxSupport = 1/n/Rawpixelsize(1);
DqySupport = 1/n/Rawpixelsize(2);
DqzSupport = 1/numz/Rawpixelsize(3);
QXSupport = ((1:n)-floor(n/2)-1)*DqxSupport;
QYSupport = ((1:n)-floor(n/2)-1)*DqySupport;
QZSupport = ((1:numz)-floor(numz/2)-1)*DqzSupport;
[qx,qy,qz] = meshgrid(QXSupport,QYSupport,QZSupport);
q0 = OTFmaskrefmed/lambda;
NAl = ReceiveNA/lambda;
NBl = sqrt(q0^2-NAl^2);
clear DqxSupport DqySupport QXSupport QYSupport QZSupport

% wide-field mask
qpar = sqrt((qx).^2+(qy).^2);
axialcutoff = sqrt(q0^2-(qpar-NAl).^2)-NBl;
axialcutoff = single(qpar<=2*NAl).*axialcutoff;
OTFmask = axialcutoff+DqzSupport/2 >= abs(qz);
OTFmask = (qpar<=2*NAl).*OTFmask;

% side-lobe mask
% pz=pixelsizeZ*numz*2*Pzrefmed*(1-cosd(60))/lambda; % OTF's area 2*r*(1-cosα)
ps=round(pixelsizeZ*numz*2*Pzrefmed*cosd(60)./lambda); % distance between sidelobe and OTF center: 2*r*cosα
r=(2*q0*360)/315;
OTFmask1=(qpar<=2*NAl).*(sqrt(qx.^2+qy.^2+qz.^2)<=r);
OTFmask1(:,:,zmedian-ps:zmedian+ps)=0;
OTFmask=OTFmask+OTFmask1; % I2M mask
% imwritestack(gather(OTFmask),'0I2M.tif')
clearvars -except OTFmask n numz xmedian

% 读取非相干照明pattern
ill=load('Functions\Illumination spectrum.mat');

% 2阶
ill_2=zeros(n,n,numz);
ill_2(xmedian,xmedian,:)=ill.A0_save(:,1);
Mask_4PS_2=convn(ill_2,OTFmask,'same');
Mask_4PS_2(Mask_4PS_2<0.1) = 0;
Mask_4PS_2(Mask_4PS_2~=0) = 1;
% imwritestack(gather(Mask_6B_2),'0ill_0.tif')

% 1阶
ill_1=zeros(n,n,numz);
ill_1(xmedian,xmedian,:)=ill.A0_save(:,2);
Mask_4PS_1=convn(ill_1,OTFmask,'same');
Mask_4PS_1(Mask_4PS_1<0.1) = 0;
Mask_4PS_1(Mask_4PS_1~=0) = 1;

% 0阶
ill_0=zeros(n,n,numz);
ill_0(xmedian,xmedian,:)=ill.A0_save(:,3);
Mask_4PS_0=convn(ill_0,OTFmask,'same');
Mask_4PS_0(Mask_4PS_0<0.1) = 0;
Mask_4PS_0(Mask_4PS_0~=0) = 1;

% 频谱拼接
Mask_4PS=zeros(n,n,numz*3);
Mask_4PS(:,:,0*numz+1:numz*1)=Mask_4PS_0; % m=0
Mask_4PS(:,:,1*numz+1:numz*2)=Mask_4PS_2; % m=2
Mask_4PS(:,:,2*numz+1:numz*3)=Mask_4PS_1; % m=1
Mask_4PS(Mask_4PS<0.1) = 0;
Mask_4PS(Mask_4PS~=0) = 1;
% imwritestack(gather(Mask_4PS),'Mask_4PS.tif')

end