clear;
lambda = 1;
mu = 1;
format = "Gauss_Seidel";
%format = "SSOR";
syms x y;
u = [-sin(pi*x)^2*sin(2*pi*y); sin(2*pi*x)*sin(pi*y)^2];
%u = [(exp(sin(pi*x))-1)*(exp(sin(pi*y))-1); sin(pi*x)*sin(pi*y)];
divu = diff(u(1), x) + diff(u(2), y);
f = mu * (diff(diff(u, x), x) + diff(diff(u, y), y)) + (lambda + mu) * [diff(divu, x); diff(divu, y)];
f = -f;
ux = matlabFunction(simplify(diff(u, x)));
uy = matlabFunction(simplify(diff(u, y)));
f = matlabFunction(simplify(f));
[phi, phix, phiy, T] = init_fespace;
error = zeros(1, 4);
rate = zeros(1, 4);
iter = zeros(1, 4);
time = zeros(1, 4);
gm = decsg([3; 4; 0; 1; 1; 0; 0; 0; 1; 1]);
for i = 1: 4
    if i == 1
        [p, e, t] = initmesh(gm, 'Hmax', 1/8);
    else
        [p, e, t] = refinemesh(gm, p, e, t);
    end
    [vertices, mesh] = init_mesh(p, t);
    not_bdr = prod(vertices .* (1 - vertices), 2);
    free = not_bdr ~= 0;
    [A, b] = FEM(vertices, mesh, lambda, mu, f, phi, phix, phiy, T);
    c = zeros(size(b));
    tic;
    [c(free), iter(i)] = eval(format+"(A(free, free), b(free), c(free), Inf)");
%    c(free) = A(free, free) \ b(free);
    time(i) = toc;
    error(i) = calc_err(vertices, mesh, lambda, mu, ux, uy, phix, phiy, T, c);
    if i > 1
        rate(i) = log2(error(i-1)/error(i));
    end
end
fprintf("error & %.3e & %.3e & %.3e & %.3e\\\\\n", error(1), error(2), error(3), error(4));
fprintf("rate & & %.2f & %.2f & %.2f\\\\\n", rate(2), rate(3), rate(4));
fprintf("iter & %d & %d & %d & %d\\\\\n", iter(1), iter(2), iter(3), iter(4));
fprintf("time & %.2f & %.2f & %.2f & %.2f\\\\\n", time(1), time(2), time(3), time(4));