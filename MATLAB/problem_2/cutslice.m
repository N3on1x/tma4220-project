n = 1000;
P = getSlice(n, 3*pi/2);
x = P(:,1);
y = P(:,2);
dt = delaunay(x,y);
eps = 3/4*pi/(2*n);

% Remove triangles 
tri2del = zeros(1, length(dt));
for i = 1:length(dt)
    k = dt(i,:);
    ce = (P(k(1),:)+P(k(2),:)+P(k(3),:))/3; % Calculate the centroid
    if (ce(1) > 0 && ce(2) < 0) % Check if it lies in the fourth quadrant
        tri2del(i) = i;
    end
end
tri2del = nonzeros(tri2del);
dt(tri2del,:) = [];

% Locate the boundry nodes
boundry = zeros(1, n);
for i = 1:n
    if ((abs(1 - norm(P(i,:))) < eps) || (abs(y(i)) < eps && x(i) >= 0) || (abs(x(i)) < eps && y(i) <= 0))
        boundry(i) = i;
    end
end
boundry = nonzeros(boundry);

%triplot(dt,x,y);
z = x.*y.*sin(2*pi*(x.^2+y.^2));
trisurf(dt,x,y,z)