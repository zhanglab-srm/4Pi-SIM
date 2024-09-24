% caculate h£¨r£©* g£¨r£©and replace 'conv2' with Fourier transform to improve computational speed
function f = dft(h, g)
	K = size(h);
	L = size(g);
	f = zeros(K+L);
	f = ifftn((fftn(h, K+L) .* fftn(g, K+L)),K+L);
	f = (f(1:K(1)+L(1)-1,1:K(2)+L(2)-1));
end