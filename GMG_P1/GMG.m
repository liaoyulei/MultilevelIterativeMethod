function c_h = GMG(vertices_h, mesh_h, l, A_h, b_h, c_h, format)
not_bdr = prod(vertices_h .* (1 - vertices_h), 2);
free = not_bdr ~= 0;
if l == 1
    c_h(free) = A_h(free, free) \ b_h(free);
else
    c_h(free) = eval(format+"(A_h(free, free), b_h(free), c_h(free), 5)");
    Nv_h = size(vertices_h, 1);
    Nt_h = size(mesh_h, 1);
    Nt_l = Nt_h / 4;
    mesh_l(:, 1) = min(mesh_h(1: Nt_l, [1, 3, 5]), [], 2);
    mesh_l(:, 2) = mesh_l(:, 1) + 1;
    mesh_l(:, 3) = min(mesh_h(Nt_l+1: 2*Nt_l, [1, 3, 5]), [], 2);
    mesh_l(:, 4) = mesh_l(:, 3) + 1;
    mesh_l(:, 5) = min(mesh_h(2*Nt_l+1: 3*Nt_l, [1, 3, 5]), [], 2);
    mesh_l(:, 6) = mesh_l(:, 5) + 1;
    Nv_l = max(mesh_l, [], 'all');
    vertices_l = vertices_h(1: Nv_l, :);
    Q = zeros(Nv_l, Nv_h);
    for k = 1: 3*Nt_l
        i = min(mesh_h(k, [1, 3, 5]));
        for j = mesh_h(k, [1, 3, 5])
            if j == i
                Q(i, j) = 1;
                Q(i+1, j+1) = 1;
            else
                Q(i, j) = 1/2;
                Q(i+1, j+1) = 1/2;
            end
        end
    end
    Q = sparse(Q);
    A_l = Q * A_h * Q';
    b_l = Q * (b_h - A_h * c_h);
    c_l = zeros(size(b_l));
    c_l = GMG(vertices_l, mesh_l, l-1, A_l, b_l, c_l, format);
    c_l = GMG(vertices_l, mesh_l, l-1, A_l, b_l, c_l, format);
    c_h = c_h + Q' * c_l;
    c_h(free) = eval(format+"(A_h(free, free), b_h(free), c_h(free), 5)");
end