% Coordinate reordering
function [tmp_zuobiaox,tmp_zuobiaoy]=zub_change(zuobiaox,zuobiaoy)
    tmp_zuobiaox = zuobiaox;
    tmp_zuobiaoy = zuobiaoy;
    Ctn = zuobiaox(1,1);
    tmp_zuobiaox([1,4,7],:) = [];
    tmp_zuobiaoy([1,4,7],:) = [];
    for chi = 2:2:size(tmp_zuobiaox,1)
        tmt = tmp_zuobiaox(chi,2);
        tmp_zuobiaox(chi,2) = tmp_zuobiaox(chi,1);
        tmp_zuobiaox(chi,1) = tmt;
        tmt = tmp_zuobiaoy(chi,2);
        tmp_zuobiaoy(chi,2) = tmp_zuobiaoy(chi,1);
        tmp_zuobiaoy(chi,1) = tmt;
    end
    tmp_zuobiaox=tmp_zuobiaox';
    tmp_zuobiaox=tmp_zuobiaox(:);
    tmp_zuobiaox=reshape([reshape(tmp_zuobiaox,2,6);[Ctn,Ctn,Ctn,Ctn,Ctn,Ctn]],1,18);
    tmp_zuobiaox(:,[6,12,18])=[];
    tmp_zuobiaoy=tmp_zuobiaoy';
    tmp_zuobiaoy=tmp_zuobiaoy(:);
    tmp_zuobiaoy=reshape([reshape(tmp_zuobiaoy,2,6);[Ctn,Ctn,Ctn,Ctn,Ctn,Ctn]],1,18);
    tmp_zuobiaoy(:,[6,12,18])=[];
end