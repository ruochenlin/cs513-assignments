function A = A3(epsilon, m)
A = epsilon * eye(m);
for i = 1 : m - 1
	for j = i + 1 : m
		A(i, j) = 1; 
	end
end

end
