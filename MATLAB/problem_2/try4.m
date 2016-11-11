close all; clc;

f = @(x,y) 16*pi^2*x.*y.*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2))...
    - 24*x.*y.*pi*cos(2*pi*(x.^2+y.^2));

n = 500;
[p, tri, edge] = getSlice(n, 3*pi/2);

A = zeros(n);
b = zeros(n,1);
[nk, np] = size(tri);

for i = 1:nk
    k = tri(i,:); % Nodes of the current element
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i]
    phi = @(x,y) ([x y 1]*coeff_mat)';
    J_phi = coeff_mat(1:2,:)'; % Jacobian of the vector function phi

    area = quadrature2D(p(k(1),:), p(k(2),:), p(k(3),:), 1, @(x,y) 1); % Area of current element.
    A_local = area*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix
    b_local = quadrature2D(p(k(1),:), p(k(2),:), p(k(3),:), 4, @(x,y) phi(x,y)*f(x,y));
    b(k') = b(k') + b_local;
end

% Boundry conditions
for i = 1:size(edge)
    node = edge(i,1);
    A(node,:) = zeros(1,n);
    A(node,node) = 1;
    b(node) = 0;
end

u_h = A\b; % Solve the system
x = p(:,1);
y = p(:,2);

figure;
trisurf(tri, x, y ,u_h)
title('FEM');

figure;
u = x.*y.*sin(2*pi*(x.^2+y.^2));
trisurf(tri, x, y, u)
title('Analytical solution');