close all; clc;

f = @(x,y) 16*pi^2*x*y*(x^2+y^2)*sin(2*pi*(x^2+y^2))...
    + 3*cos(2*pi*(x^2+y^2));

n = 10;
[p, tri, edge] = getDisk(n);

A = zeros(n);
b = zeros(n,1);
nk = 0; % Number of elements
np = 0; % Number of points per element
[nk, np] = size(tri);

J_phi = [1 0; 0 1; -1 -1]; % Jacobian of the phis
phi_dot = J_phi*(J_phi'); % Quadratic form of the Jacobian
A_local = zeros(np);
area = 0;
for i = 1:nk
    k = tri(i,:);
    p1 = p(k(1),:)';
    p2 = p(k(2),:)';
    p3 = p(k(3),:)';
    area = abs(det([(p2-p1) (p3-p1)]))/2; % Area of the current element
    A_local = area/2*phi_dot;
    A(k,k) = A(k,k) + A_local; % Build A
end