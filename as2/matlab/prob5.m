m = 2;
n=zeros(20,1);
while m <= 21
	A = randn(m,m);
	b = randn(m,1);
	[x,n(m-1)] = linearsysqr(A,b);
	m = m + 1;
end
fit(linspace(2,21,20)',n,'poly3')
fit(linspace(2,21,20)',n,'poly4')

