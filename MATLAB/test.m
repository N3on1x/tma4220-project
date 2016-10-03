close all;
%load seamount;
n = 100; % Number of mesh points
P = getDisk(n);
x = P(:,1);
y = P(:,2);
tri = delaunay(x,y);
z = x.*y.*sin(2*pi*(x.^2+y.^2));
trisurf(tri,x,y,z) % Plot the function z evaluated in the mesh points
%triplot(tri,x,y) % Plot the mesh in 2D