A = [2, -2; 1, 1];
rho = max(abs(eig(A)))
norm_1 = norm(A, 1)
norm_2 = norm(A, 2)
norm_inf = norm(A, Inf)
[U, Sigma, V] = svd(A)
