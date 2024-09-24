% 3D window
function wd = Window(nz, nx, ny, etaxy, etaz)
% windowing to avoid edge effects in FFT
padsize = 0;
x = 1:(ny+2*padsize);
y = (1:(nx+2*padsize))';
z = 1:(nz+2*padsize);
sigmaXY = etaxy;
mask = single(repmat(sigmoid(sigmaXY*(x-padsize)) - sigmoid(sigmaXY*(x-ny-padsize-1)), nx+2*padsize, 1) .* repmat(sigmoid(sigmaXY*(y-padsize)) - sigmoid(sigmaXY*(y-nx-padsize-1)), 1, ny+2*padsize));
wd = zeros(nx, ny, nz);
sigmaZ = etaz;
if nz == 1
    wd(:, :) = (mask.^3);
else
    wz = sigmoid(sigmaZ*(z-padsize)) - sigmoid(sigmaZ*(z-nz-padsize-1));
    for i = 1:nz
        wd(:, :, i) = (mask.^3) .* wz(i).^25;
    end
end
end