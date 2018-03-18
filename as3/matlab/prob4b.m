function [p, r] = prob4b(f, a, b, m, k)
% Pre: a < b; f is a function handle; m > k and are both integers
% Post: p is the coefficient vector of polynomial, r is the residue
T = linspace(a, b, m);
target = f(T)';
A = prob4a(T, k);
[Q, R] = qr(A);
target_tilde = Q' * target;
% use back substitution to solve for p
p = zeros(k + 1, 1);
for i = k + 1 : -1 : 1
	for j = k + 1 : -1 : i + 1
		target_tilde(i) = target_tilde(i) - R(i, j) * p(j);
	end
	p(i) = target_tilde(i) / R(i, i);
end

fprintf('With m=%d, k=%d, f(x) is fitted with the polynomial\n',m,k);
fprintf('%f', p(1));
for i =  1 : k
	fprintf(' + %f x^%d', p(i + 1), i);
end
r = norm(A * p - target);
fprintf('\nThe norm of residue is %f.\n', norm(A * p - target));
