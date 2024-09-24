% Use 3D OTF with pixelsizeZ numz to generate 3D OTF with pixelsizeZ_OTF numz_OTF in frequency space
function OTF_all2 = GenerateOTF_Fre_revised(OTF_all1, pixelsizeZ1, pixelsizeZ2, OTFz2, OTFxy2)
nort = 3;
[OTFx1, OTFy1, numz] = size(OTF_all1);
OTFz1 = numz/nort;
ratio = pixelsizeZ2/pixelsizeZ1;
OTF_all2 = zeros(OTFxy2, OTFxy2, OTFz2*nort);
if rem(ratio, 1) ~= 0
    warndlg('Z step must be an integer multiple of 100nm');
    return;
end
if ratio == 1  
    if OTFxy2 == OTFx1 && OTFxy2 == OTFy1 && OTFz2 == OTFz1
        OTF_all2 = OTF_all1;
    else
        for mi = 1:nort
            OTF_all2(:,:,(mi-1)*OTFz2+1:mi*OTFz2) = imresize3(OTF_all1(:,:,(mi-1)*OTFz1+1:mi*OTFz1), [OTFxy2, OTFxy2, OTFz2], 'linear');
        end
    end
else       
    MedainZ= floor(OTFz1/2)+1;
    OTFResampleZ = floor(OTFz1/ratio);
    if mod(OTFResampleZ,2) == 0
        OTFResampleZ = OTFResampleZ + 1;
    end
    OTFAllTemp = zeros(OTFx1, OTFy1, OTFResampleZ*nort);
    Z2 = MedainZ - floor(OTFResampleZ/2) : MedainZ + floor(OTFResampleZ/2);
    CycleNum = ratio - 1;
    for mi = 1:nort
        OTF = OTF_all1(:,:,(mi-1)*OTFz1+1:mi*OTFz1);        
        ResampleIndexZ = 1:OTFz1;
        ResampleOTF = OTF(:,:,ResampleIndexZ(Z2));
        for iCycleNum = 1:CycleNum
           ResampleIndexZ = mod(ResampleIndexZ-OTFResampleZ, OTFz1);
           ResampleIndexZ(ResampleIndexZ ==0 ) = OTFz1;
           IndexZ = ResampleIndexZ(Z2);
           [OTFMaxIndexZ,MaxIndexZ] = max(IndexZ);           
           if MaxIndexZ == numel(Z2)
               ResampleOTF = ResampleOTF + OTF(:,:,IndexZ);
           else
               [OTFMinIndexZ,MinIndexZ] = min(ResampleIndexZ(Z2));
               ResampleOTF(:,:,1:MaxIndexZ) = ResampleOTF(:,:,1:MaxIndexZ) + OTF(:,:,IndexZ(1):OTFMaxIndexZ);
               ResampleOTF(:,:,MinIndexZ:end) = ResampleOTF(:,:,MinIndexZ:end) + OTF(:,:,OTFMinIndexZ:IndexZ(end));
           end
        end
        OTFAllTemp(:,:,(mi-1)*OTFResampleZ+1:mi*OTFResampleZ) = ResampleOTF;
    end
    clear OTF_all1
    if ~(OTFxy2 == OTFx1 && OTFxy2 == OTFy1 && OTFz2 == OTFz1)
        for mi = 1:nort
            OTF_all2(:,:,(mi-1)*OTFz2+1:mi*OTFz2) = imresize3(OTFAllTemp(:,:,(mi-1)*OTFResampleZ+1:mi*OTFResampleZ), [OTFxy2, OTFxy2, OTFz2], 'linear');
        end
        clear OTFAllTemp
    end
end
end