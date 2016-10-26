close all; clc;

%%
sphere;

%%
% tetramesh(tet,p);
n_w = 100;
grid1  = linspace(-1,1,n_w);
[X, Y, Z] = meshgrid(grid1,grid1,grid1);
F = scatteredInterpolant(x,y,z,u_h);
V = F(X,Y,Z);
for i = 1:n_w
    for j = 1:n_w
        for k = 1:n_w
            if grid1(i)^2+grid1(j)^2+grid1(k)^2 > 1
                V(i,j,k) = nan;
            end
        end
    end
end
%%
P = patch(isosurface(X,Y,Z,V,-1/2));

%%
xlabel('x');
ylabel('y');
zlabel('z');
P.FaceColor = 'red';
P.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud