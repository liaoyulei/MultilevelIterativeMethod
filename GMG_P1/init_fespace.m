function  [phi, phix, phiy, T] = init_fespace
syms lambda1 lambda2 lambda3 x1 x2 x3 y1 y2 y3
xi1 = x2 - x3;
xi2 = x3 - x1;
xi3 = x1 - x2;
eta1 = y2 - y3;
eta2 = y3 - y1;
eta3 = y1 - y2;
T2 = det([1, x1, y1; 1, x2, y2; 1, x3, y3]);
lambda1x = eta1 / T2;
lambda2x = eta2 / T2;
lambda3x = eta3 / T2;
lambda1y = -xi1 / T2;
lambda2y = -xi2 / T2;
lambda3y = -xi3 / T2;
phi = matlabFunction([lambda1, 0; 0, lambda1; lambda2, 0; 0, lambda2; lambda3, 0; 0, lambda3]);
phix = matlabFunction([lambda1x, 0; 0, lambda1x; lambda2x, 0; 0, lambda2x; lambda3x, 0; 0, lambda3x]);
phiy = matlabFunction([lambda1y, 0; 0, lambda1y; lambda2y, 0; 0, lambda2y; lambda3y, 0; 0, lambda3y]);
T = matlabFunction(simplify(T2 / 2));
end