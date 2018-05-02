fid = fopen('prob3.dat', 'w');
B = rand(30,10);
C = rand(10);
C(1:2,:) = C(1:2,:) / 10;
for i = 0 : 5
	C(1:2, :) = C(1:2, :) * 10;
	A = B * C;
	[U, S, V] = svd(A);
	u = U(:, 1:2);
	fprintf(fid, '%2d & %4.2f & %4.5e \n', i, norm(A - u*u'*A), norm(u*u'*A));
end
