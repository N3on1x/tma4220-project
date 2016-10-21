close all; clc;

f = @(x,y) 16*pi^2*x.*y.*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2))...
    - 24*x.*y.*pi*cos(2*pi*(x.^2+y.^2));
g_x = @(x) -x.*sin(2*pi*x.^2);
g_y = @(y)  y.*sin(2*pi*y.^2);

n = 1000;
[p, tri, edge] = getSlice(n, 3*pi/2);

A = zeros(n);
b = zeros(n,1);
[nk, np] = size(tri);

% Split the boundry into the Dirichlet part and the Neumann part.
edg = edge(:,1);
ned = length(edg); % Number of edges
n_n_x = zeros(1,ned); % Nodes on the x-axis on the Neumann boundry
n_n_y = zeros(1,ned); % Nodes on the y-axis on the Neumann boundry
eps = 1e-12; % A small value. Because of rounding errors the points on the axes are not always identically zero.
for i = 1:ned
    p_edg = p(edg(i),:); % Coordinates of current edge
    % Check if a node is on the Neumann boundry
    if (abs(p_edg(2)) < eps && p_edg(1) >= 0)
        n_n_x(i) = edg(i);
    end
    if (abs(p_edg(1)) < eps && p_edg(2) <= 0)
        n_n_y(i) = edg(i);
    end
end
n_n_x = n_n_x(n_n_x ~= 0); % Remove zeros
n_n_y = n_n_y(n_n_y ~= 0);
n_d = setdiff(edg', [n_n_x n_n_y]); % Nodes on the Dirichlet boundry

for i = 1:nk
    k = tri(i,:); % Nodes of the current element
    p_k = p(k,:); % Coordinates of the nodes
    phi_sys = [p_k ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i]
    phi = @(x,y) ([x y 1]*coeff_mat)';
    J_phi = coeff_mat(1:2,:)'; % Jacobian of the vector function phi

    area = quadrature2D(p_k(1,:), p_k(2,:), p_k(3,:), 1, @(x,y) 1); % Area of current element.
    A_local = area*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix
    b_local = quadrature2D(p_k(1,:), p_k(2,:), p_k(3,:), 4, @(x,y) phi(x,y)*f(x,y));
    b(k') = b(k') + b_local;
    
    % Enforce Neumann boundry condtions
    for comb = nchoosek(1:np,np-1)' % Combinations of all the nodes in the triangle
        nodes = k(comb);
        if     sum(ismember(nodes, n_n_x)) == np-1
            % The edge is on the x-axis of the Neumann boundry
            sub_mat = zeros(np, np-1);
            sub_mat(comb,:) = eye(np-1);
            sub_phi = @(x) sub_mat'*phi(x,0);
            b(nodes) = b(nodes) + quadrature1D(p(nodes(1),1)',p(nodes(2),1)', 4, @(x) sub_phi(x)*g_x(x));
        elseif sum(ismember(nodes, n_n_y)) == np-1
            sub_mat = zeros(np, np-1);
            sub_mat(comb,:) = eye(np-1);
            sub_phi = @(y) sub_mat'*phi(0,y);
            b(nodes) = b(nodes) + quadrature1D(p(nodes(1),2)',p(nodes(2),2)', 4, @(y) sub_phi(y)*g_y(y));
        end
    end
    
end

% Enforce Dirichlet boundry conditions
for i = n_d
    A(i,:) = zeros(1,n);
    A(i,i) = 1;
    b(i) = 0;
end

u_h = A\b; % Solve the system

x = p(:,1);
y = p(:,2);

figure;
trisurf(tri, x, y ,u_h)
title('FEM');