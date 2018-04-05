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
	fprintf('Before: \n');
	U
	u_temp_row = U(k, k : end);
	U(k, k : end) = U(max_index(1), k : end);
	U(max_index(1), k : end) = u_temp_row;
	u_temp_col = U(:, k);
	U(:, k) = U(:, max_index(2));
	U(:, max_index(2)) = u_temp_col;
	fprintf('After: \n')
	U
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
