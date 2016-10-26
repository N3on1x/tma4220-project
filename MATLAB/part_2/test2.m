
close all; clc; clear all;
n = 100;
p = linspace(0, 1, n)';
tet = [(1:n-1)' (2:n)'];
edge = [1;n];
[n, nd] = size(p);
[nk, np] = size(tet);

A = zeros(n);
M = zeros(n);

for i = 1:nk
    k = tet(i,:); % Nodes of the current element
    phi_sys = [p(k,:) ones(np, 1)];
    coeff_mat = phi_sys\eye(np); % Collumn i is [a_i; b_i]
    phi = @(x,y,z) ([x 1]*coeff_mat)'; % phi vector
    J_phi = coeff_mat(1:nd,:)'; % Jacobian of the vector function phi

    volume = quadrature1D(p(k(1),:), p(k(2),:), 1, @(x) 1); % Volume of current element.
    A_local = volume*J_phi*(J_phi'); % Local stiffness matrix
    A(k,k) = A(k,k) + A_local; % Build the stiffness matrix
    
    M_local = quadrature1D(p(k(1),:), p(k(2),:), 4, @(x) phi(x)*(phi(x)'));
    M(k,k) = M(k,k) + M_local;
end

edg = unique(edge); % Boundry nodes
bnd_d = 225; % Constant Dirichlet boundry condition
intn = setdiff(1:n,edg); % Internal nodes

c = A(intn,edg)*bnd_d*ones(length(edg),1);
A(edg,:) = [];
A(:,edg) = [];
M(edg,:) = [];
M(:,edg) = [];

%%

u_h = 20*ones(length(intn),1); % Initial condition
h = 1e-3; % Time step
t_end = 10; % End time
for i = 0:h:t_end
    plot(p(intn),u_h);
    waitforbuttonpress
    %u_h = (M+h/2*A)\((M-h/2*A)*u_h-h*c); % CN
    u_h = (M+h*A)\(M*u_h-h*c);
end
plot(p(intn), u_h);