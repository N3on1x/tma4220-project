%%
% This script changes the surface id number in the fifth column of the tri
% data to simplify.

clear all, close all, clc;

data_edge = dlmread('../../Meshes/heatsink055_tri.m');

[nrow,ncol] = size(data_edge);

robinEdgeIndexes = find(ismember(data_edge(:,5), 6));
otherEdgeIndexes = setdiff(1:nrow,robinEdgeIndexes);

data_edge(robinEdgeIndexes,5) = 1;
data_edge(otherEdgeIndexes,5) = 2;

dlmwrite('../../Meshes/heatsink055_tri_new.m',data_edge,' ');

%%
data_edge = dlmread('../../Meshes/heatsink105_tri.m');

[nrow,ncol] = size(data_edge);

robinEdgeIndexes = find(ismember(data_edge(:,5), 6));
otherEdgeIndexes = setdiff(1:nrow,robinEdgeIndexes);

data_edge(robinEdgeIndexes,5) = 1;
data_edge(otherEdgeIndexes,5) = 2;

dlmwrite('../../Meshes/heatsink105_tri_new.m',data_edge,' ');