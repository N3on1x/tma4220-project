clc; close all; clear all;
% Load the mesh
data_nodes = dlmread('../../Meshes/heatsink105_nodes.m');
data_tet = dlmread('../../Meshes/heatsink105_tetr.m');
data_edge = dlmread('../../Meshes/heatsink105_tri_new.m');
p = data_nodes(:,2:end);
tet = data_tet(:,1:4);
edge = data_edge(:,1:3);
[n, nd] = size(p); % Number of points n in dimension nd
[nk, np] = size(tet); % Number of elements nk and points per element np

% Physical properties of the problem
h_tr = 10; % Transfer coefficient between copper and air
k_cu = 401; % Thermal conductivity of copper.
T_amb = 20; % Air temperature
bnd_d = 80; % Constant Dirichlet boundry condition

A = zeros(n);
b = zeros(n,1);

for i = 1:nk
    k = tet(i,:); % Nodes of the current element
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    J_phi = coeff_mat(1:nd,:)'; % Jacobian of the vector function phi

    volume = quadrature3D(p(k(1),:), p(k(2),:), p(k(3),:), p(k(4),:), 1, @(x,y,z) 1); % Volume of current element.
    A_local = k_cu/h_tr*volume*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix 
end

% Split the edge into the different boundries
edge_d_ind   = find(data_edge(:,5) == 1);
edge_r_ind_z = find(data_edge(:,5) == 2);
edge_r_ind_x = find(data_edge(:,5) == 4);
edge_r_ind_y = find(data_edge(:,5) == 3);

edge_d   = edge(edge_d_ind,:);
edge_r_x = edge(edge_r_ind_x,:);
edge_r_y = edge(edge_r_ind_y,:);
edge_r_z = edge(edge_r_ind_z,:);

% Enfore Robin boundry condition
for i = 1:length(edge_r_z)
    k = edge_r_z(i,:); % Current triangle on the Robin boundry
    phi_sys = [[p(k,:); zeros(1,nd)] ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    coeff_mat(:,np) = [];
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    pt = p(k,:); % Coordinates
    nv = cross(pt(2,:)-pt(1,:),pt(3,:)-pt(1,:)); % Normal vector of the plane
    z = @(x,y) (dot(nv,pt(1,:))-nv(1)*x-nv(2)*y)/nv(3); % Equation of a plane
    A_local = quadrature2D(pt(1,[1 2]), pt(2,[1 2]), pt(3,[1 2]), 4, @(x,y) phi(x,y,z(x,y))*(phi(x,y,z(x,y))')*norm(nv)/abs(nv(3)));
    A(k,k) = A(k,k) + A_local;
    b_local = T_amb*quadrature2D(pt(1,[1 2]), pt(2,[1 2]), pt(3,[1 2]), 4, @(x,y) phi(x,y,z(x,y))*norm(nv)/abs(nv(3)));
    b(k) = b(k) + b_local;
end

for i = 1:length(edge_r_y)
    k = edge_r_y(i,:); % Current triangle on the Robin boundry
    phi_sys = [[p(k,:); -ones(1,nd)] ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    coeff_mat(:,np) = [];
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    pt = p(k,:); % Coordinates
    nv = cross(pt(2,:)-pt(1,:),pt(3,:)-pt(1,:)); % Normal vector of the plane
    y = @(x,z) (dot(nv,pt(1,:))-nv(1)*x-nv(3)*z)/nv(2); % Equation of a plane
    A_local = quadrature2D(pt(1,[1 3]), pt(2,[1 3]), pt(3,[1 3]), 4, @(x,z) phi(x,y(x,z),z)*(phi(x,y(x,z),z)')*norm(nv)/abs(nv(2)));
    A(k,k) = A(k,k) + A_local;
    b_local = T_amb*quadrature2D(pt(1,[1 3]), pt(2,[1 3]), pt(3,[1 3]), 4, @(x,z) phi(x,y(x,z),z)*norm(nv)/abs(nv(2)));
    b(k) = b(k) + b_local;
end

for i = 1:length(edge_r_x)
    k = edge_r_x(i,:); % Current triangle on the Robin boundry
    phi_sys = [[p(k,:); -ones(1,nd)] ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    coeff_mat(:,np) = [];
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    pt = p(k,:); % Coordinates
    nv = cross(pt(2,:)-pt(1,:),pt(3,:)-pt(1,:)); % Normal vector of the plane
    x = @(y,z) (dot(nv,pt(1,:))-nv(2)*y-nv(3)*z)/nv(1); % Equation of a plane
    A_local = quadrature2D(pt(1,[2 3]), pt(2,[2 3]), pt(3,[2 3]), 4, @(y,z) phi(x(y,z),y,z)*(phi(x(y,z),y,z)')*norm(nv)/abs(nv(1)));
    A(k,k) = A(k,k) + A_local;
    b_local = T_amb*quadrature2D(pt(1,[2 3]), pt(2,[2 3]), pt(3,[2 3]), 4, @(y,z) phi(x(y,z),y,z)*norm(nv)/abs(nv(1)));
    b(k) = b(k) + b_local;
end

% Enforce Dirichlet boundry conditions
d_nodes = unique(edge_d); % Dirichlet boundry nodes
not_d_nodes = setdiff(1:n,d_nodes); % The rest of the nodes
c = A(not_d_nodes,d_nodes)*bnd_d*ones(length(d_nodes),1);
A(d_nodes,:) = [];
A(:,d_nodes) = [];
b(d_nodes) = [];

u_h = bnd_d*ones(n,1);
u_h(not_d_nodes) = A\(b-c);