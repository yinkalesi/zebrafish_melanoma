%% findIntersect_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/19/19
%  Project: Fish Image Analysis
%  find the first pixel in the line between two points on the path
function [inter,dist] = findIntersect_v1_0(p1,p2,map)

[bounds,line] = getBounds_v1_0(p1,p2,map);
dist = find(bounds,1,'first');
inter = line(dist,:);