close all;clc;clear all;
data_nodes = dlmread('data/cake_nodes.dat');
data_tet = dlmread('data/cake_tetr.dat');
data_edge = dlmread('data/cake_tri.dat');

p = data_nodes(:,2:4);
edge = data_edge(:,1:3);
tet = data_tet(:,1:4);
%es = find(ismember(data_edge(:,5), [122,88,105]));
tr = triangulation(edge,p);
trimesh(tr)