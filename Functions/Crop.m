%调节矩阵大小，从f中裁剪出和g维数一致的矩阵（裁剪四周，保留中间）
function f = Crop(f, g)
	K = size(g);
	N = size(f);
	L = ceil((N-K) / 2);
	v = colonvec(L+1, L+K);
	f = f(v{:});
end