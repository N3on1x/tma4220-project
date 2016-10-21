% Analytical solution for problem 2 on the unit disk.
%close all;
n = 1000; % Number of mesh points
P = getSlice(n,3*pi/2);
x = P(:,1);
y = P(:,2);
tri = delaunay(x,y);

z = x.*y.*sin(2*pi*(x.^2+y.^2));
trisurf(tri,x,y,z) % Plot the function z evaluated in the mesh points
%triplot(tri,x,y) % Plot the mesh in 2D