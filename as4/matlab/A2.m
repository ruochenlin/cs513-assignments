function A = A2(m)
for i = 1 : m
	for j = 1 : m
		A(i, j) = i^(j - 1);
	end
end
