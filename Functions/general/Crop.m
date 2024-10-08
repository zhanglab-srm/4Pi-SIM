% crop images
function f = Crop(f, g)
	K = size(g);
	N = size(f);
	L = ceil((N-K) / 2);
	v = colonvec(L+1, L+K);
	f = f(v{:});
end