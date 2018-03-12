function result = prob1(m, n, s)
% Pre: m, n are positive integers, s is a 2*2 matrix with integer entries
[X1, X2] = meshgrid(1:1:m, 1:1:n);
candidates = [X1(:)'; X2(:)'];
x = s \ candidates;
result = candidates(:, all(x == round(x)));
