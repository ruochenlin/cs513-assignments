function [Q, S] = Q1(A, b, n)
	m = length(A);
	Q(:, 1) = b / norm(b);
	for i = 1 : n
		% Assuming we already have Q(:, 1:i)
		v = A * Q(:, i);
		for j = 1 : i
			s(j) = Q(:, j)' * v;
			v = v - s(j) * Q(:, j);
		end
		norm_v = norm(v);
		if norm_v < 1e-14
			S = S(1:end-1, :);
			break;
		end
		S(:, i) = s; 
		if i ~= n
			S(i+1, i) = norm_v;
			Q(:, i+1) = v / norm_v;
		end
	end
	
end
