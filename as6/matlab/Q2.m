function err = Q2(A, b, n)
	warning off;
	m = length(A);
	x_matlab = A \ b;
	Q(:, 1) = b / norm(b);
	for i = 1 : n
		% Assuming we already have Q(:, 1:i)
		v = A * Q(:, i);
		clear s;
		for j = 1 : i
			s(j) = Q(:, j)' * v;
			v = v - s(j) * Q(:, j);
		end
		norm_v = norm(v);
		% If we have found optimal x in K_i
		if norm_v < 1e-14
			break;
		end
		S(:, i) = s;
		S(i+1, i) = norm_v;
		if i ~= n
			Q(:, i+1) = v / norm_v;
		end
	end
	b_tilde = zeros(length(S), 1);
	b_tilde(1) = norm(b);
	y = S \ b_tilde;
	x_gmres = Q * y;
	err = norm(x_matlab - x_gmres);
end
