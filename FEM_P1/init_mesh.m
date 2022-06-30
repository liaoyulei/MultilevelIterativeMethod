function [vertices2, mesh2] = init_mesh(p, t)
Nv = size(p, 2);
Nt = size(t, 2);
vertices = p';
mesh = t(1: 3, :)';
N = 2 * Nv;
vertices2 = zeros(N, 2);
vertices2(1: 2: N, :) = vertices;
vertices2(2: 2: N, :) = vertices;
mesh2 = zeros(Nt, 6);
mesh2(:, 1: 2: 6) = 2 * mesh - 1;
mesh2(:, 2: 2: 6) = 2 * mesh;
end