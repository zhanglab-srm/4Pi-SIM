%调节矩阵大小函数
function v = colonvec(m, M)
	n = numel(m);
	N = numel(M);
	K = max(n, N);
	v = cell(K, 1);
	if n == 1
		m = m * ones(K, 1);
	elseif N == 1
		M = M * ones(K, 1);
	end
	for k = 1:K
		v{k} = m(k):M(k);
	end	
end