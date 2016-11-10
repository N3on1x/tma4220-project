close all;
[x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
v = x.*exp(-x.^2-y.^2-z.^2);
xslice = NaN; 
yslice = 0.5; 
zslice = NaN;
slice(x,y,z,v,xslice,yslice,zslice)
xlabel('x');
ylabel('y');
zlabel('z');