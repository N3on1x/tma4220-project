%% 
%%% Load the mesh
data_nodes = dlmread('data/cake_nodes.dat');
data_tet = dlmread('data/cake_tetr.dat');
data_edge = dlmread('data/cake_tri.dat');
p = data_nodes(:,2:end);
tet = data_tet(:,1:4);
tetr_material = data_tet(:,5);
edge = data_edge(:,1:3);
[n, nd] = size(p);
[nk, np] = size(tet);

%tetramesh(tet,p)
%% 

A = zeros(n);
M = zeros(n);

for i = 1:nk
    k = tet(i,:); % Nodes of the current element
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i; c_i; d_i]
    phi = @(x,y,z) ([x y z 1]*coeff_mat)'; % phi vector
    J_phi = coeff_mat(1:nd,:)'; % Jacobian of the vector function phi

    volume = quadrature3D(p(k(1),:), p(k(2),:), p(k(3),:), p(k(4),:), 1, @(x,y,z) 1); % Volume of current element.
    A_local = volume*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix
    
    M_local = quadrature3D(p(k(1),:), p(k(2),:), p(k(3),:), p(k(4),:), 5, @(x,y,z) phi(x,y,z)*(phi(x,y,z)'));
    M(k,k) = M(k,k) + M_local;
end