function [apod_func] = apodization(nx, nz, plong, fc,pixelsizeZ,Exwave)
    RI = 1.33;
    Emwave = 525;
    
    pz = pixelsizeZ*nz*RI/Exwave;
    OTF_1z = pixelsizeZ*nz*RI/Emwave;
    fz = (pz + OTF_1z)*2*0.85; % 4PS apodization
    
    if mod(nz,2) == 1 % z layers number is odd
        [xv, yv, zv] = meshgrid(-(2*nx)/2:(2*nx)/2-1, -(2*nx)/2:(2*nx)/2-1, -(floor(nz/2)):floor(nz/2));
    else
        [xv, yv, zv] = meshgrid(-(2*nx)/2:(2*nx)/2-1, -(2*nx)/2:(2*nx)/2-1, -(nz)/2:(nz)/2-1);
    end
    
    k_rxy = sqrt(xv.^2+yv.^2);
    k_rz = sqrt(zv.^2);
    k_max_xy = plong(1)+fc;
    apod_func = cos(pi*sqrt(k_rxy.^2./k_max_xy.^2+0*k_rz.^2./fz.^2)/2); % lateral apodization + GUASSIAN
    apod_func = apod_func.*cos(pi*sqrt(0*k_rxy.^2./k_max_xy.^2+k_rz.^2./fz.^2)/2); % axial apodization + GUASSIAN
    indi =  k_rxy > k_max_xy ;
    apod_func(indi) = 0;
    indii =  k_rz > fz ;
    apod_func(indii) = 0;
end
