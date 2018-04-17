% In this script, we compare the performance of two algorithms
% that gives the number of eigenvalues in (0,1) of a symmetric
% tridiagonal matrix. To avoid the overhead of indexing large 
% arrays, we always try to save the values of matrix elements 
% that can be reused.
clear all;
BASE = 100;
count = 15;
for i = 1 : count
	m = i * BASE;
	t_lu_i = 0;
	t_det_i = 0;
	for j = 1 : 1000
		% construct a random symmetric tridiagonal matrix
		x = rand(m, 1);
		y = rand(m-1, 1);
		A1 = diag(x) + diag(y, 1) + diag(y, -1);
		A2 = A1;
		B1 = A1 - eye(m);
		B2 = B1;
		
		% TESTING
		% q = eig(A1);
		% n = length(q(q > 0 & q < 1));

		% find the number of negative eigenvalues of A1
		% with lu factorization of tridiagonal matrix
		clear a_kp1_k;clear a_kp1_kp1;clear a_k_k; n_lu_0 = 0;
		tic;
		a_k_k = A1(1,1);
		if a_k_k < 0
			n_lu_0 = 1;
		end
		for k = 1 : m - 1
			a_kp1_k = A1(k+1,k);
			a_kp1_kp1 = A1(k+1, k+1);
			a_kp1_kp1 = a_kp1_kp1 - a_kp1_k^2 / a_k_k;
			a_k_k = a_kp1_kp1;
			if a_kp1_kp1 < 0
				n_lu_0 = n_lu_0 + 1;
			end
		end		
		% find the number of negative eigenvalues of B1
		% with lu factorization of tridiagonal matrix
		clear b_kp1_k;clear b_kp1_kp1;clear b_k_k; n_lu_1 = 0;
		b_k_k = B1(1,1);
		if b_k_k < 0
			n_lu_1 = 1;
		end
		for k = 1 : m - 1
			b_kp1_k = B1(k+1,k);
			b_kp1_kp1 = B1(k+1, k+1);
			b_kp1_kp1 = b_kp1_kp1 - b_kp1_k^2 / b_k_k;
			b_k_k = b_kp1_kp1;
			if b_kp1_kp1 < 0
				n_lu_1 = n_lu_1 + 1;
			end
		end
		n_lu = n_lu_1 - n_lu_0;
		t_lu_i = t_lu_i + toc;

		% find the number of negative eigenvalues of A1
		% by calculating determinants of principal minors
		% of a tridiagonal matrix
		clear d_km1; clear d_km2; clear d_k; n_det_0 = 0;
		tic;
		d_km2 = 1;
		d_km1 = A2(1,1);
		if d_km1 < 0
			n_det_0 = 1;
		end
		for k = 2 : m
			d_k = A2(k,k) * d_km1 - A2(k, k-1)^2 * d_km2;
			% if sign changes, it means there's a negative eigenvalue
			if xor(d_k < 0, d_km1 < 0)
				n_det_0 = n_det_0 + 1;
			end
			d_km2 = d_km1;
			d_km1 = d_k;
		end
		% find the number of negative eigenvalues of B1
		% by calculating determinants of principal minors
		% of a tridiagonal matrix
		clear d_km1; clear d_km2; clear d_k; n_det_1 = 0;
		d_km2 = 1;
		d_km1 = B2(1,1);
		if d_km1 < 0
			n_det_1 = 1;
		end
		for k = 2 : m
			d_k = B2(k,k) * d_km1 - B2(k, k-1)^2 * d_km2;
			% if sign changes, it means there's a negative eigenvalue
			if xor(d_k < 0, d_km1 < 0)
				n_det_1 = n_det_1 + 1;
			end
			d_km2 = d_km1;
			d_km1 = d_k;
		end
		n_det = n_det_1 - n_det_0;
		t_det_i = t_det_i + toc;
		% if n_det ~= n_lu | n ~= n_det
		% 	fprintf('!!!!!!  ');
		% end
		% fprintf('%4d: %4d, %4d, %4d\n', i*BASE, n, n_lu, n_det);
	end
	t_lu(i) = t_lu_i;
	t_det(i) = t_det_i;
end
plot(linspace(BASE, count*BASE, count), t_lu, linspace(BASE, count*BASE, count), t_det);
xlabel('Size of Matrix');
ylabel('Time (s)');
legend('LU', 'Determinant');
