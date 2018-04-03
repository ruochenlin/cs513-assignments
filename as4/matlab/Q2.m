function [C, e1, e2, e3] = Q2(A, x)
C = cond(A);
b = A * x;
m = length(A);

% Solve Ax = b with original LU
if issymmetric(A)
	% reuse the symm_lu routine if A is symmetric
	[L, U] = symm_lu(A);
else
	U1 = A; L1 = eye(m);
	for k = 1 : m - 1
		L1(k + 1 : end, k) = U1(k + 1 : end, k) / U1(k, k);
		U1(k + 1 : end, k) = 0;
		for i = k + 1 : m
			U1(i, k + 1 : end) = U1(i, k + 1 : end) - L1(i, k) * U1(k, k + 1 : end);
		end
	end
end
% Use backsubstitution to solve L1*y1 = b1 and then U1*x1 = y1 
b1 = b;
y1 = zeros(m, 1);
for i = 1 : m
	for j = 1 : i - 1
		b1(i) = b1(i) - L1(i, j) * y1(j);
	end
	y1(i) = b1(i) / L1(i, i);
end
x1 = zeros(m, 1);
for i = m : -1 : 1
	for j = m : -1 : i + 1
		y1(i) = y1(i) - U(i, j) * x1(j);
	end
	x1(i) = y1(i) / U(i, i);
end
e1 = x - x1;

% Use Matlab lu routine
[L2, U2, P] = lu(A);
b2 = P * b;
% Use backsubstitution to solve L2*y2 = b2 and then U2*x2 = y2
y2 = zeros(m, 1);
for i = 1 : m
    for j = 1 : i - 1 
        b2(i) = b2(i) - L2(i, j) * y2(j);
    end
    y2(i) = b2(i) / L2(i, i);
end
x2 = zeros(m, 1);
for i = m : -1 : 1
    for j = m : -1 : i + 1
        y2(i) = y2(i) - U(i, j) * x2(j);
    end
    x2(i) = y2(i) / U(i, i);
end
e2 = x - x2;

% Use LU factorization with full pivoting