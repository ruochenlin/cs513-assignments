function [C, e1, e2, e3] = Q2(A, x)
C = cond(A);
b = A * x;
m = length(A);

% Solve Ax = b with original LU
if issymmetric(A)
	% use the symm_lu routine if A is symmetric
	[L1, U1] = symm_lu(A);
else
	[L1, U1] = lu_no_pivot(A);
end
% Solve L1*y1 = b1 and then U1*x1 = y1 
b1 = b;
y1 = zeros(m, 1);
for i = 1 : m
	for j = 1 : i - 1
		b1(i) = b1(i) - L1(i, j) * y1(j);
	end
	y1(i) = b1(i);
end
x1 = zeros(m, 1);
for i = m : -1 : 1
	for j = m : -1 : i + 1
		y1(i) = y1(i) - U1(i, j) * x1(j);
	end
	x1(i) = y1(i) / U1(i, i);
end
e1 = norm(x - x1);

% Use Matlab lu routine
[L2, U2, P2] = lu(A);
b2 = P2 * b;
% Solve L2*y2 = b2 and then U2*x2 = y2
y2 = zeros(m, 1);
for i = 1 : m
    for j = 1 : i - 1 
        b2(i) = b2(i) - L2(i, j) * y2(j);
    end
    y2(i) = b2(i);
end
x2 = zeros(m, 1);
for i = m : -1 : 1
    for j = m : -1 : i + 1
        y2(i) = y2(i) - U2(i, j) * x2(j);
    end
    x2(i) = y2(i) / U2(i, i);
end
e2 = norm(x - x2);

% Use LU factorization with full pivoting
[P3, Q, L3, U3] = lu_full_pivot(A);
b3 = P3 * b;
% Solve L3 * z = b3, then U3 * y3 = z, finally x3 = Q * y3 
z = zeros(m, 1);
for i = 1 : m
    for j = 1 : i - 1
        b3(i) = b3(i) - L3(i, j) * z(j);
    end
    z(i) = b3(i);
end
y3 = zeros(m, 1);
for i = m : -1 : 1
    for j = m : -1 : i + 1
        z(i) = z(i) - U3(i, j) * y3(j);
    end
    y3(i) = z(i) / U3(i, i);
end
x3 = Q * y3;
e3 = norm(x - x3);

end

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

function [P, Q, L, U] = lu_full_pivot(A)
% Pre: A is a non-singular square matrix
% Post: P*A*Q = L*U
SMALL_NUM = 1e-15;
m = length(A);
U = A; L = eye(m); P = eye(m); Q = eye(m);

for k = 1 : m - 1
	% find the entry with max abs value in subblock
	max_index = [k, k];
	max_entry = abs(U(k, k));
	for i = k : m
		for j = k : m
			if i == k && j == k
				continue;
			end
			if max_entry < abs(U(i, j))
				max_entry = abs(U(i, j));
				max_index = [i, j];
			end
		end
	end
	
	% Permute P, Q, L, U according to the position of pivot
	if max_index(1) ~= k
		P([k, max_index(1)], :) = P([max_index(1), k], :);
	end
	% p_temp = P(k, :);
	% P(k, :) = P(max_index(1), :);
	% P(max_index(1), :) = p_temp;
	q_temp = Q(:, k);
	Q(:, k) = Q(:, max_index(2));
	Q(:, max_index(2)) = q_temp; 
	if k > 1
		l_temp = L(k, 1 : k - 1);
		L(k, 1 : k - 1) = L(max_index(1), 1 : k - 1);
		L(max_index(1), 1 : k - 1) = l_temp;
	end
	u_temp_row = U(k, k : end);
	U(k, k : end) = U(max_index(1), k : end);
	U(max_index(1), k : end) = u_temp_row;
	u_temp_col = U(:, k);
	U(:, k) = U(:, max_index(2));
	U(:, max_index(2)) = u_temp_col;

	% check if dividing by 0
	if abs(U(k, k)) < SMALL_NUM
		err = MException('flag:DivideByZero', 'Denominator is (almost) zero!');
		U
		throw(err);
	end 
	
	% Update L, U
	L(k + 1 : end, k) = U(k + 1 : end, k) / U(k, k);
	U(k + 1 : m, k) = 0;
	for i = k + 1 : m
		U(i, k + 1 : end) = U(i, k + 1 : end) - L(i, k) * U(k, k + 1 : end);
	end
end

end
