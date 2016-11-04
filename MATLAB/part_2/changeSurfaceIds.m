%%
% This script changes the surface id number in the fifth column of the tri
% data to simplify.

clear all, close all, clc;

data_edge = dlmread('../../Meshes/heatsink055_tri.m');

[nrow,ncol] = size(data_edge);

dirchletEdgeIndexes = find(ismember(data_edge(:,5), 6));
xySurfaceIndexes = find(ismember(data_edge(:,5),[134,132,130,128,118,96,74]));
xzSurfaceIndexes = find(ismember(data_edge(:,5),[124,122,117,109,95,87,73,65]));
yzSurfaceIndexes = find(ismember(data_edge(:,5),[126,120,113,105,91,83,69,61]));

data_edge(dirchletEdgeIndexes,5) = 1;
data_edge(xySurfaceIndexes,5) = 2;
data_edge(xzSurfaceIndexes,5) = 3;
data_edge(yzSurfaceIndexes,5) = 4;

dlmwrite('../../Meshes/heatsink055_tri_new.m',data_edge,' ');

%%
data_edge = dlmread('../../Meshes/heatsink105_tri.m');

[nrow,ncol] = size(data_edge);

dirchletEdgeIndexes = find(ismember(data_edge(:,5), 6));
xySurfaceIndexes = find(ismember(data_edge(:,5),[134,132,130,128,118,96,74]));
xzSurfaceIndexes = find(ismember(data_edge(:,5),[124,122,117,109,95,87,73,65]));
yzSurfaceIndexes = find(ismember(data_edge(:,5),[126,120,113,108,91,83,69,61]));

data_edge(dirchletEdgeIndexes,5) = 1;
data_edge(xySurfaceIndexes,5) = 2;
data_edge(xzSurfaceIndexes,5) = 3;
data_edge(yzSurfaceIndexes,5) = 4;

dlmwrite('../../Meshes/heatsink105_tri_new.m',data_edge,' ');

%%
data_edge = dlmread('../../Meshes/heatsink_h405_gmf005_tri.m');

[nrow,ncol] = size(data_edge);

dirchletEdgeIndexes = find(ismember(data_edge(:,5), 6));
xySurfaceIndexes = find(ismember(data_edge(:,5),[134,132,130,128,118,96,74]));
xzSurfaceIndexes = find(ismember(data_edge(:,5),[124,122,117,109,95,87,73,65]));
yzSurfaceIndexes = find(ismember(data_edge(:,5),[126,120,113,108,91,83,69,61]));

data_edge(dirchletEdgeIndexes,5) = 1;
data_edge(xySurfaceIndexes,5) = 2;
data_edge(xzSurfaceIndexes,5) = 3;
data_edge(yzSurfaceIndexes,5) = 4;

dlmwrite('../../Meshes/heatsink_h405_gmf005_tri_new.m',data_edge,' ');