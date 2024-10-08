function [phase, contrast] = Fitting_histogram(band0,moveOTF,temp_band,OTF) 
    % This function performs histogram statistics on the overlapped bands to
    % extract phase and contrast.
    
    addpath('Functions/package_emd/EMDs')
    precision1 = 0.02;
    precision2 = 0.02;
    f = -pi:precision1:pi;
    ff = 0:precision2:0.6;
    overlapRegion = band0.*moveOTF;
    cm_con = temp_band./(overlapRegion+eps);
    cm_ang = temp_band./(overlapRegion+eps);
    angcm = angle(cm_ang);  % initial phase
    abscm = abs(cm_con);    % contrast
    replc6_ang = moveOTF.*OTF;
    replc6_con = replc6_ang;
    a_ang = replc6_ang(:,:,ceil((size(OTF,3)+1)/2));
    a_con = replc6_con(:,:,ceil((size(OTF,3)+1)/2));
    b_ang = a_ang(:);
    b_con = a_con(:);
    b_ang(b_ang~=0) = 1;
    b_con(b_con~=0) = 1;
    c = angcm(:,:,ceil((size(OTF,3)+1)/2));
    cc = abscm(:,:,ceil((size(OTF,3)+1)/2));
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
    g = histc(e,f);
    gg = histc(ee,ff);
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
    % emd 
    imf=emd(g);
    imf2=emd(gg);
    if max(g)<50     
        g=sum(imf(5:end,:),1);
    else
        g=sum(imf(4:end,:),1);
    end
    gg=sum(imf2(1:end,:),1);
    h=find(g==max(g));
    if h > num-Index+1               
        phase = -pi+precision1*(mean(h)-(num-Index+1));
    else
        phase = -pi+precision1*(mean(h)+Index-1);
    end
    contrast = mean(ee(ee<1)); 
end