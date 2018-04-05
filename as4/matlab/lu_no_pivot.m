function [L, U] = lu_no_pivot(A)
% Pre: A is a nonsingular square matrix
% Post: A = L*U
SMALL_NUM = 1e-15;
if issymmetric(A)
	% reuse the symm_lu routine if A is symmetric
	[L, U] = symm_lu(A);
else
	m = length(A);
	U = A; L = eye(m);
	for k = 1 : m - 1
		% check if dividing by 0
		if abs(U(k, k)) < SMALL_NUM
            err = MException('flag:DivideByZero', 'Denominator is (almost) zero!');
            throw(err);
		end
		% Update L, U
		L(k + 1 : end, k) = U(k + 1 : end, k) / U(k, k);
		U(k + 1 : end, k) = 0;
		for i = k + 1 : m
			U(i, k + 1 : end) = U(i, k + 1 : end) - L(i, k) * U(k, k + 1 : end);
		end
	end
end

end
