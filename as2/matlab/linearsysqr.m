function [x, n] = linearsysqr(A,b)
% a function that solves linear systems with QR factorization
% and counts the number of operations
m = max(size(A));
n = 0;
k = 1;
while k < m
	y = zeros(m,1);
	if k > 1
		y(1:k-1) = A(1:k-1,k); n = n + k - 1;
	end
	y(k) = norm(A(k:end, k)); n = n + 1;
	w = A(:, k) - y; n = n + m;
	w = w/norm(w); n = n + 3*m - 1;
	b = b - 2*(w'*b)*w; n = n + 4*m;
	A (:, k) = y; n = n + m;
	i = k + 1;
	while i < m + 1
		A (:, i) = A(:, i) - 2 * (w'*A(:,i))*w; n = n + 4 * m;
		i = i + 1;
	end
	k = k + 1;
end
x = zeros(m, 1); n = n + m;
k = m;
while k > 0
	i = m;
	while i > k
		b(k) = b(k) - A(k,i)*x(i); n = n + 2;
		i = i - 1;
	end
	x(k) = b(k) / A(k,k); n = n + 1;
	k = k - 1;
end
