low = 3; high = 30;
t = zeros(high - low + 1, 1);
for m = low : high
	for i = 1 : 100
		A = rand(m,m);
		A = A' + A;
		tic;
		[L,U] = symm_lu(A);
		t(m - low + 1) = t(m - low + 1) + toc;
	end
end
t = t / 100;
