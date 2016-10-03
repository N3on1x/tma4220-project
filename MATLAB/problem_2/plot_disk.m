% Plots the the mesh on the unit disk
n = 21; % Number of nodes
P = getDisk(n);
tri = delaunayTriangulation(P);
%triplot(tri);

el = size(tri);
el = el(1); % Number of elements
A = zeros(n);

for i = 1:el
    k = tri(i,:);
    p1 = P(k(1),:);
    p2 = P(k(2),:);
    p3 = P(k(3),:);
end
