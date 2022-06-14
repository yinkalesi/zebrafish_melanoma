%% findPathDistance_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/20/19
%  Project: Fish Image Analysis
%  find (approximate) shortest distance along path between two pixels

function [dist] = findPathDistance_v1_0(p1,p2,map)
    anchors = findAnchorPoints_v1_0(p1,p2,map);
    dz = [anchors;p2]-[p1;anchors];
    ds = sqrt(sum(dz.^2,2));
    dist = sum(ds);
end