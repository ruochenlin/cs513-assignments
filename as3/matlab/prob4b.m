function p = prob4b(f, a, b, m, k)
% Pre: a < b; f is a function handle; m > k and are both integers
T = linspace(a, b, m);
target = f(T)';
A = prob4a(T, k);
[Q, R] = qr(A);
target_tilde = Q' * target;
% using '\' to solve a linear system of equations, not least squares 
p(1 : k + 1, 1) = R(1 : k + 1, 1 : k + 1) \ (target_tilde(1 : k + 1));
fprintf('f(x) is fitted with the polynomial\n');
fprintf('%f', p(1));
for i =  1 : k
	fprintf(' + %f x^%d',p(i+1), i);
end
fprintf('\nThe norm of error is %f.\n', norm(A * p - target));
