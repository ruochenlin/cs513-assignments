function A = A1(m)
for i = 1 : m
	for j = 1 : m
		A(i, j) = 1 / (i + j - 1); 
	end
end
