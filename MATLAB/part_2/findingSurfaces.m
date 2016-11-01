clc, clear all, close all;

data_edge = dlmread('../../Meshes/heatsink055_tri.m');
data_nodes = dlmread('../../Meshes/heatsink055_nodes.m');

p = data_nodes(:,2:end);
edge = data_edge(:,1:3);

surfaceIds = unique(data_edge(:,5));

for i = 1:length(surfaceIds)
    id = surfaceIds(i);
    surfEdges = find(ismember(data_edge(:,5), id));
    figure();
    trimesh(triangulation(edge(surfEdges,:),p));
    zlim([-0.5,1.5]);
    xlim([0 1]);
    ylim([0 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(num2str(id));
end
