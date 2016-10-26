% Heat equation on a rod of length L with initial conditions
% f = @(x) 6*sin(pi*x/L); and homogenous dirichlet boundry conditions.
% Known analytic solution u = @(x,t) 6*sin(pi*x/L).*exp(-alpha*(pi/L)^2*t)

close all; clc;
L = 1;
alpha = 1;
n = 100;
x = linspace(0, L, n);
f = @(x) 6*sin(pi*x/L);
u = @(x,t) 6*sin(pi*x/L).*exp(-alpha*(pi/L)^2*t);

p = x';
tet = [(1:n-1)' (2:n)'];
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
A = alpha*A;

% Impose Dirichlet conditions
A([1 n],:) = [];
A(:, [1 n]) = [];
M([1 n],:) = [];
M(:, [1 n]) = [];

h = 1e-3;
x = x(2:n-1);
u_h = f(x)';
t_end = 1; % End time
for i = 0:h:t_end
    u_h = (2*M+h*A)\((2*M-h*A)*u_h); % CN
end
plot(x,u_h)
hold on;
plot(x, u(x,t_end))
legend('fem','anal')