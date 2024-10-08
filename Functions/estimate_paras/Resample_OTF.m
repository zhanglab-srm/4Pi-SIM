function OTF_new = Resample_OTF(OTF_all_raw, pixelsizeZ_raw, pixelsizeZ_new, OTFz_new, OTFxy_new)
    % This function will resample OTF according to SIM data.
    
    nort = 3;
    [OTFx_raw, OTFy_raw, numz] = size(OTF_all_raw);
    OTFz_raw = numz/nort;
    ratio = pixelsizeZ_new/pixelsizeZ_raw;
    OTF_new = zeros(OTFxy_new, OTFxy_new, OTFz_new*nort);
    if rem(ratio, 1) ~= 0
        warndlg('Z step must be an integer multiple of 40nm');
        return;
    end
    if ratio == 1  
        if OTFxy_new == OTFx_raw && OTFxy_new == OTFy_raw && OTFz_new == OTFz_raw
            OTF_new = OTF_all_raw;
        else
            for mi = 1:nort
                OTF_new(:,:,(mi-1)*OTFz_new+1:mi*OTFz_new) = imresize3(OTF_all_raw(:,:,(mi-1)*OTFz_raw+1:mi*OTFz_raw), [OTFxy_new, OTFxy_new, OTFz_new], 'linear');
            end
        end
    else       
        MedianZ= floor(OTFz_raw/2)+1;
        OTFResampleZ = floor(OTFz_raw/ratio);
        if mod(OTFResampleZ,2) == 0
            OTFResampleZ = OTFResampleZ + 1;
        end
        OTFAllTemp = zeros(OTFx_raw, OTFy_raw, OTFResampleZ*nort);
        Z2 = MedianZ - floor(OTFResampleZ/2) : MedianZ + floor(OTFResampleZ/2);
        CycleNum = ratio - 1;
        for mi = 1:nort
            OTF = OTF_all_raw(:,:,(mi-1)*OTFz_raw+1:mi*OTFz_raw);        
            ResampleIndexZ = 1:OTFz_raw;
            ResampleOTF = OTF(:,:,ResampleIndexZ(Z2));
            for iCycleNum = 1:CycleNum
               ResampleIndexZ = mod(ResampleIndexZ-OTFResampleZ, OTFz_raw);
               ResampleIndexZ(ResampleIndexZ ==0 ) = OTFz_raw;
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
        clear OTF_all_raw
        if ~(OTFxy_new == OTFx_raw && OTFxy_new == OTFy_raw && OTFz_new == OTFz_raw)
            for mi = 1:nort
                OTF_new(:,:,(mi-1)*OTFz_new+1:mi*OTFz_new) = imresize3(OTFAllTemp(:,:,(mi-1)*OTFResampleZ+1:mi*OTFResampleZ), [OTFxy_new, OTFxy_new, OTFz_new], 'linear');
            end
            clear OTFAllTemp
        end
    end
    
    OTFz_new = size(OTF_new,3)/nort;
    OTF_new(:,:,1:OTFz_new) = OTF_new(:,:,1:OTFz_new)./max(max(max(abs(OTF_new(:,:,1:OTFz_new))))); % m = 0
    OTF_new(:,:,1*OTFz_new+1:2*OTFz_new) = OTF_new(:,:,1*OTFz_new+1:2*OTFz_new)./max(max(max(abs(OTF_new(:,:,1*OTFz_new+1:2*OTFz_new))))); % m = 2
    OTF_new(:,:,2*OTFz_new+1:3*OTFz_new) = OTF_new(:,:,2*OTFz_new+1:3*OTFz_new)./max(max(max(abs(OTF_new(:,:,2*OTFz_new+1:3*OTFz_new))))); % m = 1


end