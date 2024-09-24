% 3D bhs z guassian mask
function [bhs] = bhs_4PS_theory(nx, nz, plong, fc,pixelsizeZ,wavelengh,I2M_para)
n=1.33;
NA_ex= 1; %1.15
NA_min= 1; %1.1
EmissionWavelength=525;
% I2M_para = 5;

% d = n-sqrt(n.^2-NA_ex.^2);
% d = n;
pz = pixelsizeZ*nz*n/wavelengh; % 4PS - 0
OTF_1z = pixelsizeZ*nz*n/EmissionWavelength;

fz = (pz + OTF_1z)*2*0.85; % 4PS apodization

if mod(nz,2) == 1 % z layers number is odd
    [xv, yv, zv] = meshgrid(-(2*nx)/2:(2*nx)/2-1, -(2*nx)/2:(2*nx)/2-1, -(floor(nz/2)):floor(nz/2));
else
    [xv, yv, zv] = meshgrid(-(2*nx)/2:(2*nx)/2-1, -(2*nx)/2:(2*nx)/2-1, -(nz)/2:(nz)/2-1);
end

k_rxy = sqrt(xv.^2+yv.^2);
k_rz = sqrt(zv.^2);
k_max_xy = plong(1)+fc;
bhs = cos(pi*sqrt(k_rxy.^2./k_max_xy.^2+0*k_rz.^2./fz.^2)/2); % lateral apodization +GUASSIAN
bhs = bhs.*cos(pi*sqrt(0*k_rxy.^2./k_max_xy.^2+k_rz.^2./fz.^2)/2); % axial apodization+GUASSIAN
indi = find( k_rxy > k_max_xy );
bhs(indi) = 0;
indii = find( k_rz > fz );
bhs(indii) = 0;
% imwritestack(bhs,'bhs_4PS-0.85.tif');
end
