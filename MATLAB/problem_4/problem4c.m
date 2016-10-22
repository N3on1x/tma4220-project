close all; clc;

sphere2;
% tetramesh(tet,p);
n_w = 100;
grid  = linspace(-1.1,1.1,n_w);
[X, Y, Z] = meshgrid(grid,grid,grid);
V = sin(2*pi*X).*sin(2*pi*Y).*sin(2*pi*Z);
P = patch(isosurface(X,Y,Z,V,1/2));
isonormals(X,Y,Z,V,P)
P.FaceColor = 'red';
P.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud