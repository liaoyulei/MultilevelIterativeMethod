function error = calc_err(vertices, mesh, lambda, mu, ux, uy, phix, phiy, T, c)
ld = [1/3, 1/3, 1/3];
error = 0;
normu = 0;
Nt = size(mesh, 1);
a = @(ux, uy, vx, vy) lambda * (ux(1) + uy(2)) * (vx(1) + vy(2)) + 2 * mu * (ux(1) * vx(1) + uy(2) * vy(2)) + mu * (uy(1) + ux(2)) * (vy(1) + vx(2));
for k = 1: Nt
    v1 = vertices(mesh(k, 1), :);
    v2 = vertices(mesh(k, 3), :);
    v3 = vertices(mesh(k, 5), :);
    vx = phix(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    vy = phiy(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    x = v1(1) * ld(1) + v2(1) * ld(2) + v3(1) * ld(3);
    y = v1(2) * ld(1) + v2(2) * ld(2) + v3(2) * ld(3);
    ex = ux(x, y) - vx' * c(mesh(k, :));
    ey = uy(x, y) - vy' * c(mesh(k, :));
    error = error + jacobi * a(ex, ey, ex, ey);
    normu = normu + jacobi * a(ux(x, y), uy(x, y), ux(x, y), uy(x, y));
end
error = (error / normu)^(1/2);
end