A = [1, 2, 3, 4, 5, 6;
2, -2, -1, 1, 3, 5;
3, -1, 3, 5, 7, 9;
4, 1, 5, -4, 2, 8;
5, 3, 7, 2, 5, -10;
6, 5, 9, 8, -10, -6]
[VA, LambdaA] = eig(A)
[VA2, LambdaA2] = eig(A^2)
