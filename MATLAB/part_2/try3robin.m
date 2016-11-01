%%
clc; clear all;
% Load the mesh
data_nodes = dlmread('data/cake_nodes.dat');
data_tet = dlmread('data/cake_tetr.dat');
data_edge = dlmread('data/cake_tri.dat');
p = data_nodes(:,2:end);
tet = data_tet(:,1:4);
edge = data_edge(:,1:3);
[n, nd] = size(p);
[nk, np] = size(tet);

%tetramesh(tet,p)
%% 

A = zeros(n);
b = zeros(n,1);
h = 1;
k_cu = 1;
T_amb = 20; % Air temperature
bnd_d = 80; % Constant Dirichlet boundry condition

for i = 1:nk
    k = tet(i,:); % Nodes of the current element
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    J_phi = coeff_mat(1:nd,:)'; % Jacobian of the vector function phi

    volume = quadrature3D(p(k(1),:), p(k(2),:), p(k(3),:), p(k(4),:), 1, @(x,y,z) 1); % Volume of current element.
    A_local = k_cu/h*volume*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix 
end

% Robin boundry conditions
edg_ind_r = find(ismember(data_edge(:,5), [42,59,27,76,96,113,130]));
edg_r = edge(edg_ind_r,:);


for i = 1:length(edg_r)
    k = edg_r(i,:);
    phi_sys = [[p(k,:); zeros(1,np-1)] ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    coeff_mat(:,np) = [];
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    pt = p(k,:); % Coordinates
    nv = cross(pt(2,:)-pt(1,:),pt(3,:)-pt(1,:)); % Normal vector of the plane
    z = @(x,y) (dot(nv,pt(1,:))-nv(1)*x-nv(2)*y)/nv(3); % Equation of a plane
    A_local = quadrature2D(pt(1,1:2), pt(2,1:2), pt(3,1:2), 4, @(x,y) phi(x,y,z(x,y))*(phi(x,y,z(x,y))')*norm(nv)/abs(nv(3)));
    A(k,k) = A(k,k) + A_local;
    b_local = T_amb*quadrature2D(pt(1,1:2), pt(2,1:2), pt(3,1:2), 4, @(x,y) phi(x,y,z(x,y))*norm(nv)/abs(nv(3)));
    b(k) = b(k) + b_local;
end

% Dirichlet boundry conditions
edg_ind_d = setdiff(1:length(edge), edg_ind_r);
edg_d = edge(edg_ind_d,:);

edg = unique(edg_d); % Boundry nodes
not_d = setdiff(1:n,edg); % Elements that are not on the Dirichlet boundry
c = A(not_d,edg)*bnd_d*ones(length(edg),1);
A(edg,:) = [];
A(:,edg) = [];
b(edg) = [];

u_h = bnd_d*ones(n,1);
u_h(not_d) = A\(b-c);