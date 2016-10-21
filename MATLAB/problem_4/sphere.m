% Gir riktig svar, sikker p? hvorfor.

close all; clc;

f = @(x,y,z) 12*pi^2*sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);

n = 500;
[p, tet, edge] = getSphere(n);
[n, n_d] = size(p);

A = zeros(n);
b = zeros(n,1);
[nk, np] = size(tet);

for i = 1:nk
    k = tet(i,:); % Nodes of the current element
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    phi = @(x,y,z) ([x y z 1]*coeff_mat)';
    J_phi = coeff_mat(1:n_d,:)'; % Jacobian of the vector function phi

    volume = quadrature3D(p(k(1),:), p(k(2),:), p(k(3),:), p(k(4),:), 1, @(x,y,z) 1); % Volume of current element.
    A_local = volume*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix
    b_local = quadrature3D(p(k(1),:), p(k(2),:), p(k(3),:), p(k(4),:), 4, @(x,y,z) phi(x,y,z).*f(x,y,z));
    b(k') = b(k') + b_local;
end

x = p(:,1);
y = p(:,2);
z = p(:,3);
u = sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);

% Dirichlet boundry conditions
edg = unique(edge);
n_edg = length(edg);
A(edg,:) = zeros(n_edg,n);
A(edg,edg) = eye(n_edg);
b(edg) = u(edg);

u_h = A\b; % Solve the system

err = norm(u-u_h, inf);