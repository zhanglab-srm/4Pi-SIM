function [tmp_coodinatex,tmp_coodinatey] = Transform_coodinate(coodinatex,coodinatey)
    % This function transforms coodinate of wave vectors.
    tmp_coodinatex = coodinatex;
    tmp_coodinatey = coodinatey;
    Ctn = coodinatex(1,1);
    tmp_coodinatex([1,4,7],:) = [];
    tmp_coodinatey([1,4,7],:) = [];
    for chi = 2:2:size(tmp_coodinatex,1)
        tmt = tmp_coodinatex(chi,2);
        tmp_coodinatex(chi,2) = tmp_coodinatex(chi,1);
        tmp_coodinatex(chi,1) = tmt;
        tmt = tmp_coodinatey(chi,2);
        tmp_coodinatey(chi,2) = tmp_coodinatey(chi,1);
        tmp_coodinatey(chi,1) = tmt;
    end
    tmp_coodinatex = tmp_coodinatex';
    tmp_coodinatex = tmp_coodinatex(:);
    tmp_coodinatex = reshape([reshape(tmp_coodinatex,2,6);[Ctn,Ctn,Ctn,Ctn,Ctn,Ctn]],1,18);
    tmp_coodinatex(:,[6,12,18])=[];
    tmp_coodinatey = tmp_coodinatey';
    tmp_coodinatey = tmp_coodinatey(:);
    tmp_coodinatey = reshape([reshape(tmp_coodinatey,2,6);[Ctn,Ctn,Ctn,Ctn,Ctn,Ctn]],1,18);
    tmp_coodinatey(:,[6,12,18])=[];
end