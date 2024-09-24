function [angle6, c6, histogram] = SIM_3D_angle_EMD(spzhongxin,replch,youhuatest,H,orti,Pindex) 
addpath('Functions/package_emd/EMDs')
jdjd = 0.02;
jdc = 0.02;
f=-pi:jdjd:pi;
ff=0:jdc:0.6;
chongdie=spzhongxin.*replch;
cm_con = youhuatest./(chongdie+eps);
cm_ang = youhuatest./(chongdie+eps);
angcm = angle(cm_ang);
abscm = abs(cm_con);
replc6_ang = replch.*H;
replc6_con = replc6_ang;
a_ang = replc6_ang(:,:,ceil((size(H,3)+1)/2));
a_con = replc6_con(:,:,ceil((size(H,3)+1)/2));
b_ang = a_ang(:);
b_con = a_con(:);
b_ang(b_ang~=0) = 1;
b_con(b_con~=0) = 1;
c = angcm(:,:,ceil((size(H,3)+1)/2));
cc = abscm(:,:,ceil((size(H,3)+1)/2));
% imwritestack(abs(abscm),"abscm.tif")
d = c(:);
dd = cc(:);
b_ang(:,2) = d;
b_con(:,2) = dd;
b_ang((b_ang(:,1)==0),:) = [];
b_con((b_con(:,1)==0),:) = [];
e = b_ang(:,2);
ee = b_con(:,2);
e=e';
ee=ee';
testee=ee.*ee(end:-1:1);
testee=testee.^0.5;
ee=testee;
g=histc(e,f);
gg=histc(ee,ff);
% expanding
num = numel(g);
gyantou = [g(1:end-1) g(end-1) g(1:end-1) g(end-1)]; 
if sum(g) == 0                                       
    Index=1;
else
    [~, Index] = min(g(1:end-1));                   
    gyantou = gyantou(Index:Index+num-1);
    clear g
    g = gyantou;
end
histogram = g;
% figure;plot(g);
% hold on
imf=emd(g);
imf2=emd(gg);
if max(g)<50     
    g=sum(imf(5:end,:),1);
else
    g=sum(imf(4:end,:),1);
end
gg=sum(imf2(1:end,:),1);
% plot(g,'r');
% title('The Histogram of Parameter Angle')
% hold off;
% figure
% plot(gg);
% title('The Histogram of Parameter C')
h=find(g==max(g));
hh=(find(gg==max(gg)));
if h > num-Index+1               
    angle6=-pi+jdjd*(mean(h)-(num-Index+1));
else
    angle6=-pi+jdjd*(mean(h)+Index-1);
end
% angle6=-pi+mean(h)*jdjd;
% c6=jdc*mean(hh); %EMD
c6=mean(ee(ee<1)); %average value