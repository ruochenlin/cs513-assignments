function [L,U] = symm_lu(A)
% Pre: A is a real symmetric matrix
% Post: L, U are the resulting matrices of LU-factorization of A
m = length(A);
L = eye(m);
for i = 1 : m - 1
	for j = i + 1 : m
		L(j, i) = A(j, i) / A(i, i);
		A(j, i) = 0;
		A(j, i + 1 : j) = A(j, i + 1 : j) - L(j, i) * A(i, i + 1 : j);
		if i < j - 1 
			A(i + 1 : j - 1, j) = A(j, i + 1 : j - 1)';
		end
	end
end
U = A;
