function [res] = addByColumn(vals,factor)
%% addByColumn
%  do a column-wise addition; it is straightforward to do this in Matlab
%  version 2017 but doesn't work for 2015.
%  dist: MxN matrix
%  factor: 1xN matrix

res = zeros(size(vals));
for i = 1:size(vals,2)
    res(:,i) = vals(:,i)+factor(i);
end
end