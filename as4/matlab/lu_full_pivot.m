function [P, Q, L, U] = lu_full_pivot(A)
m = length(A);
U = A; L = eye(m); P = eye(m); Q = eye(m);
for k = 1 : m - 1
	max_index = [k, k];
	max_entry = abs(U(k, k));
	for i = k + 1 : m
		for j = k + 1 : m
			if max_entry < abs(U(i, j))
				max_entry = abs(U(i, j));
				max_index = [i, j];
			end
		end
	end
	p_temp = P(k, :);
	P(k, :) = P(max_index(1), :);
	P(max_index(1), :) = p_temp;
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
	L(k + 1 : end, k) = U(k + 1 : end, k) / U(k, k);
	U(k + 1 : m, k) = 0;
	for i = k + 1 : m
		U(i, k + 1 : end) = U(i, k + 1 : end) - L(i, k) * U(k, k + 1 : end);
	end
end
