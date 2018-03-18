A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
b = [1; 0; 1];
[Q, R] = qr(A)
b0 = Q'*b
b0_tilde = b0(1:2)
A1 = R(1:2, 1:2)
A2 = R(1:2, 1:2:3)
A3 = R(1:2, 2:3)
x1_tilde = A1 \ b0_tilde;
x2_tilde = A2 \ b0_tilde;
x3_tilde = A3 \ b0_tilde;
x1 = [x1_tilde; 0]
x2 = [x2_tilde(1); 0; x2_tilde(2)]
x3 = [0; x3_tilde]
fprintf('norm(A*x1-b): \n');
norm(A*x1-b)
fprintf('norm(A*x2-b): \n');
norm(A*x2-b)
fprintf('norm(A*x3-b): \n');
norm(A*x3-b)
