function [L,U] = symm_lu(A)
% Pre: A is a real symmetric matrix
% Post: L, U are the resulting matrices of LU-factorization of A
SMALL_NUM = 1e-15;
m = length(A);
L = eye(m); U = A;
for i = 1 : m - 1
	% Check if dividing by zero
	if abs(U(i, i)) < SMALL_NUM
		err = MException('flag:DivideByZero', 'Denominator is (almost) zero!');
		throw(err);
	end
	% Update L, U using only upper triangal part of U
	L(i + 1 : end, i) = U(i, i + 1 : end)' / U(i, i);
	for j = i + 1 : m
		U(j, j : end) = U(j, j : end) - U(i, j : end) * L(j, i);
	end
	U(i + 1 : end, i) = 0;
end

end
