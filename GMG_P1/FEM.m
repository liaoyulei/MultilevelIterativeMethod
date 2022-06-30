function [A, b] = FEM(vertices, mesh, lambda, mu, f, phi, phix, phiy, T)
ld = [1/3, 1/3, 1/3];
Nv = size(vertices, 1);
Nt = size(mesh, 1);
a = @(ux, uy, vx, vy) lambda * (ux(1) + uy(2)) * (vx(1) + vy(2)) + 2 * mu * (ux(1) * vx(1) + uy(2) * vy(2)) + mu * (uy(1) + ux(2)) * (vy(1) + vx(2));
x = zeros(Nt*6^2, 1);
y = zeros(Nt*6^2, 1);
v = zeros(Nt*6^2, 1);
idx = @(k, i) (k - 1) * 6^2 + (i - 1) * 6 + (1: 6);
for k = 1: Nt
    K = zeros(6);
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    vx = phix(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    vy = phiy(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    for j = 1: 6
        for i = 1: 6
            K(i, j) = jacobi * a(vx(j, :), vy(j, :), vx(i, :), vy(i, :));
        end
    end
    for i = 1: 6
        x(idx(k, i)) = mesh(k, i);
        y(idx(k, i)) = mesh(k, :);
        v(idx(k, i)) = K(i, :);
    end
end
A = sparse(x, y, v, Nv, Nv);
b = zeros(Nv, 1);
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    v = phi(ld(1), ld(2), ld(3));
    x = v1(1) * ld(1) + v2(1) * ld(2) + v3(1) * ld(3);
    y = v1(2) * ld(1) + v2(2) * ld(2) + v3(2) * ld(3);
    b(mesh(k, :)) = b(mesh(k, :)) + jacobi * v * f(x, y);
end
end