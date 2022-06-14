function [res] = multiplyByColumn(vals,factor)
%% getNorm
%  do a column-wise multiplication; it is straightforward to do this in Matlab
%  version 2017 with '.*' operator but that doesn't work for 2015. So I
%  wrote this little guy
%  dist: MxN matrix
%  factor: 1xN matrix

res = zeros(size(vals));
for i = 1:size(vals,2)
    res(:,i) = vals(:,i)/factor(i);
end
end