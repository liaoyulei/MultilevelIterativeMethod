function [x, iter] = SSOR (A, b, x, iter_max)
omega = 1;
iter = 0;
eps = 1e-8; %||b-Ax||_2<=eps
N = length(b);
D = spdiags(diag(A), 0, N, N);
L = tril(D-A);
U = triu(D-A);
M1 = D - omega .* L;
M2 = D - omega .* U;
N1 = omega .* U + (1 - omega) .* D;
N2 = omega .* L + (1 - omega) .* D;
r = b - A * x;
while r' * r > eps^2 && iter < iter_max
    x = M1 \ (N1 * x + omega .* b);
    x = M2 \ (N2 * x + omega .* b);
    r = b - A * x;
    iter = iter + 1;
end
end