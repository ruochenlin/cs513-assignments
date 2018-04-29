clear all;
center = [4, 5, 10]';
fid = fopen('prob2.dat', 'w');
for i = 1 : 3 : 10
	fprintf(fid, '%3d  ', i);
	epsilon = 10^(-i);
	D = diag(center(randi([1,3], 100, 1)) + epsilon * (rand(100,1) - 0.5));
	P = rand(100, 100);
	A = P * D * inv(P);
	n = [5, 10, 20, 90];
	% n = 4:2:90;
	b = rand(100,1);
	for i = 1 : length(n)
		err(i) = Q2(A, b, n(i));
		fprintf(fid, '& %4.2f ', log10(err(i)));
	end
	% plot(n, log10(err));
	% hold on;
	fprintf(fid, '\\\\\n');
end
fclose(fid);
% hold off;
% legend('1', '2', '3', '4', '5', '6', '7', '8');
