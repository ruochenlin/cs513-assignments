function A = A4(epsilon, m)
A = zeros(m,m);
for i = 1 : m
	A(i, i) = epsilon;
	for j = 1 : m
		if i ~= j
			A(i, j) = 1;
		end
	end
end

end
