%% loadDistribution
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 6/24/15
%  Project: Tumor Growth, Discrete Form
%  [dist,x,t] = loadDistribution(loc)

function [dist,x,t] = loadDistribution(loc)
% loadDistribution obtains data stored in a distribution data file
% loc: string; file directory and name - should have M+1 rows and N+1
% columns
% dist: MxN matrix; contians data
% x: 1xM vector; value corresponding to rows of matrix
% t: 1xN vector; value corresonding to columns of matrix

% store data in matrix
mat = dlmread(loc);
dist = mat(2:end,2:end);
x = mat(2:end,1);
t = mat(1,2:end);



