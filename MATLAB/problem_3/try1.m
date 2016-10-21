close all; clc;

f = @(x,y) 16*pi^2*x.*y.*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2))...
    - 24*x.*y.*pi*cos(2*pi*(x.^2+y.^2));

n = 50;
[p, tri, edge] = getSlice(n, 3*pi/2);

A = zeros(n);
b = zeros(n,1);
[nk, np] = size(tri);

% Split the boundry into the Dirichlet part and the Neumann part.
ned = length(edge); % Number of edges
ind_edg_n = zeros(1,ned);
eps = 1e-12;
for i = 1:ned
    p_ed1 = p(edge(i,1),:); % Coordinates of current edge
    p_ed2 = p(edge(i,2),:);
    % Check if a node is on the Neumann boundry
    if ((abs(p_ed1(2)) < eps && p_ed1(1) >= 0) || (abs(p_ed1(1)) < eps && p_ed1(2) <= 0))...
    && ((abs(p_ed2(2)) < eps && p_ed2(1) >= 0) || (abs(p_ed2(1)) < eps && p_ed2(2) <= 0))
        ind_edg_n(i) = i;
    end
end
ind_edg_n = ind_edg_n(ind_edg_n ~= 0); % Remove zeros
ind_edg_d = 1:ned;
ind_edg_d(ind_edg_n) = [];
edge_n = edge(ind_edg_n,:);
edge_d = edge(ind_edg_d,:);


for i = 1:nk
    k = tri(i,:); % Nodes of the current element
    p_k = p(k,:); % Coordinates of the nodes
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i]
    phi = @(x,y) ([x y 1]*coeff_mat)';
    J_phi = coeff_mat(1:2,:)'; % Jacobian of the vector function phi

    area = quadrature2D(p_k(1), p_k(2), p_k(3), 1, @(x,y) 1); % Area of current element.
    A_local = area*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix
    b_local = quadrature2D(p_k(1), p_k(2), p_k(3), 4, @(x,y) phi(x,y)*f(x,y));
    b(k') = b(k') + b_local;
    
    % Enforce Neumann boundry condtions
    comb = nchoosek(1:np,np-1);
end

% Boundry conditions

% Dirichlet
for i = 1:size(edge_d)
    node = edge_d(i,1);
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