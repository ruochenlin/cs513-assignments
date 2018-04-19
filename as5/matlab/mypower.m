clear all;
load as5_q3.mat;
epsilon = 1e-2;
m = length(A);
% Find the eigenvalue with the largest amplitude
v = rand(m, 1);
v = v / norm(v);
Av = A * v;
norm_Av = norm(Av);
lambda = v' * Av;
while norm(Av - lambda * v) / norm_Av >= epsilon	
	v = Av / norm_Av;
	Av = A * v;
	norm_Av = norm(Av);
	lambda = v' * Av;
end

% Deflate (lambda,v) from A
e1 = zeros(m, 1);
e1(1) = 1;
w = v - e1;
w = w / norm(w);
wTA = w' * A;
Aw = A * w;
B = A - 2 * w * wTA - 2 * Aw * w' + 4 * (wTA * w) * w * w';
% B1 is the deflated (m-1) * (m-1) matrix
B1 = B(2:end, 2:end);
n = m - 1;
v1 = rand(n, 1);
v1 = v1 / norm(v1);
B1v1 = B1 * v1;
norm_B1v1 = norm(B1v1);
lambda_1 = v1' * B1v1;
while norm(B1v1 - lambda_1 * v1) / norm_B1v1 >= epsilon
	v1 = B1v1 / norm_B1v1;
	B1v1 = B1 * v1;
	norm_B1v1 = norm(B1v1);
	lambda_1 = v1' * B1v1;
end
