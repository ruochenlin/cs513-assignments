function A = prob4a(T, k)
% Pre: T is a row vector, and k is a positive integer
m = length(T);
% Generate A without using loop
A = ( T' * ones(1, k + 1) ) .^ ( ones(m, 1) * linspace(0, k, k + 1) );
