SMALL_NUM = 1e-15;
t = zeros(10, 2);
for b = 10 : 10 : 100
	 for a = 1 : 100
	 	A = rand(b);
	 	A = A' + A;
	 	tic;
	 	m = length(A);
	 	L = eye(m); U = A;
	 	for i = 1 : m - 1
	 		% Check if dividing by zero
	 		if abs(U(i, i)) < SMALL_NUM
	 			err = MException('flag:DivideByZero', 'Denominator is (almost) zero!');
	 			throw(err);
	 		end
	 		% Update L, U using only upper triangal part of U
	 		% L(i + 1 : end, i) = U(i, i + 1 : end)' / U(i, i);
	 		for j = i + 1 : m
				L(j, i) = U(i, j) / U(i, i);
	 			% U(j, j : end) = U(j, j : end) - U(i, j : end) * L(j, i);
				for k = j : m
					U(j, k) = U(j, k) - U(i, k) * L(j, i);
				end
				U(j, i) = 0;
	 		end
	 		% U(i + 1 : end, i) = 0;
	 	end
	 	t(b / 10, 1) = t(b / 10, 1) + toc;
	 	tic;
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
				L(i, k) = U(i, k) / U(k, k);
				U(i, k) = 0; 
				for j = k + 1 : m
					U(i, j) = U(i, j) - L(i, k) * U(k, j);
				end
	 			% U(i, k + 1 : end) = U(i, k + 1 : end) - L(i, k) * U(k, k + 1 : end);
	 		end
	 	end
	 	t(b / 10, 2) = t(b / 10, 2) + toc;
	 end
end
