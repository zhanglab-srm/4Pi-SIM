%DIRSTRETCH variable directional stretching of a complex vector
function y=dirstretch(x,dir,coef)
if nargin == 3 
  if ~isreal(dir) && all(abs(dir) ~= 1)
    warning('non unit modulus of direction argument ignored')
  end
  if isreal(dir) && any(abs(dir) ~= 1)
    dir = exp(1i*dir);
  end
  if any(dir==0)
    error('invalid zero direction argument')
  end
  tmp = abs(dir);
  tmp(tmp==0) = 1; % to avoid NaNs
  dir = dir./tmp;
end
if nargin == 2
    coef = abs(dir);
    tmp = coef;
    tmp(coef==0) = 1; % to avoid NaNs
    dir = 1./tmp.*dir;
end
rotx = x.*conj(dir);
y = coef.*dir.*real(rotx)+1i*dir.*imag(rotx);