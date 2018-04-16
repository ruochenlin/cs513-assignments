function n = eig_lu(A)
% Pre: A is a non-singular symmetric tridiagonal matrix
% Post: the number of negative eigenvalues of A returned
if A(1, 1) < 0
	n = 1;
else 
	n = 0;
end

for i = 1 : length(A) - 1
	A(i+1, i+1) = A(i+1, i+1) - A(i+1, i)^2 / A(i, i);
	if A(i+1, i+1) < 0
		n = n + 1;
	end
end

end
