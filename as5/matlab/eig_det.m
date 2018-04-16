function n = eig_det(A)
% Pre: A is nonsingular tridiagonal matrix
% Post: n is the number of negative eigenvalues of A
d(1) = 1;
d(2) = A(1,1);
if d(2) < 0
	n = 1;
else 
	n = 0;
end
for i = 2 : length(A)
	d(3) = A(i,i) * d(2) - A(i, i-1)^2 * d(1);
	if xor(d(3) < 0, d(2) < 0)
		n = n + 1;
	end
	d(1) = d(2);
	d(2) = d(3);
end

end
