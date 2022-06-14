function [coor] = getTransformPointsAsArray(tform,pos)
%% getTranformPointsAsArray
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 5/24/17
%  Project: Image Alignment

[u,v] = transformPointsInverse(tform,pos(1),pos(2));
coor = [u v];
end