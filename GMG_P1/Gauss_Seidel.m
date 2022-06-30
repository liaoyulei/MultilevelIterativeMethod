function [x, iter] = Gauss_Seidel(A, b, x, iter_max)
iter = 0;
eps = 1e-8; %||b-Ax||_2<=eps
N = length(b);
D = spdiags(diag(A), 0, N, N);
L = tril(D-A);
U = triu(D-A);
r = b - A * x;
while r' * r > eps^2 && iter < iter_max 
    x = (D - L) \ (U * x + b);
    r = b - A * x;
    iter = iter + 1;
end
end