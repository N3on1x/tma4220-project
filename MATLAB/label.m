close all;
n = 50;
P = getSlice(n, 3*pi/2);
x = P(:,1);
y = P(:,2);
dt = delaunayTriangulation(x,y);
triplot(dt);
%
% Display the Vertex and Triangle labels on the plot
hold on
vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:n)');
Hpl = text(x, y, vxlabels, 'FontWeight', 'bold', 'HorizontalAlignment',...
   'center', 'BackgroundColor', 'none');
ic = incenter(dt);
numtri = size(dt,1);
trilabels = arrayfun(@(x) {sprintf('T%d', x)}, (1:numtri)');
Htl = text(ic(:,1), ic(:,2), trilabels, 'FontWeight', 'bold', ...
   'HorizontalAlignment', 'center', 'Color', 'blue');
hold off