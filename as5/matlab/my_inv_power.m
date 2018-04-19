clear all;
load as5_q3.mat
epsilon = 1e-2;
mu = 0;
m = length(A);
% Use inverse power iteration to find the eigenvalue 
% of A with the smallest amplitude (closest to 0)
A1 = A;
if mu ~= 0
	for i = 1 : m
		A1(i, i) = A1(i, i) - mu;
	end
end
v = rand(m, 1);
v = v / norm(v);
[Q, R] = qr(A1);
y = Q \ v;
u = R \ y;
% zeta = 1 / (lambda - mu)
zeta = v' * u;
while norm(u-zeta*v) / norm(u) >= epsilon
	v = u / norm(u);
	y = Q \ v;
	u = R \ y;
	zeta = v' * u;
end
lambda = 1 / zeta + mu;

% Deflate (lambda, v) from A
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
if mu ~= 0
	for i = 1 : n
		B1(i, i) = B1(i, i) - mu;
	end
end
% Carry our inverse iteration on B1
v1 = rand(n, 1);
v1 = v1 / norm(v1);
[Q1, R1] = qr(B1);
y1 = Q1 \ v1;
u1 = R1 \ y1;
% zeta = 1 / (lambda - mu)
zeta_1 = v1' * u1;
while norm(u1 - zeta_1 * v1) / norm(u1) >= epsilon
	v1 = u1 / norm(u1);
	y1 = Q1 \ v1;
	u1 = R1 \ y1;
	zeta_1 = v1' * u1;
end
lambda_1 = 1 / zeta_1 + mu;
