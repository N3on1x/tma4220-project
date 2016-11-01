%close all;clc;clear all;
close all;
data_nodes = dlmread('data/cake_nodes.dat');
data_edge = dlmread('data/cake_tri.dat');


p = data_nodes(:,2:4);
edge = data_edge(:,1:3);
edge_bowl = edge(find(ismember(data_edge(:,5), [122,88,105])),:);
edge_top = edge(find(ismember(data_edge(:,5), [42,59,27,76,96,113,130])),:);
%es = find(ismember(data_edge(:,5), [122,88,105]));
tr = triangulation(edge_top,p);
trimesh(tr)