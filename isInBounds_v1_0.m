%% isInBounds_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/19/19
%  Project: Fish Image Analysis
%  see if pixels in line of pixels in between two points all lie on path
function [allin] = isInBounds_v1_0(p1,p2,map)

allin = all(getBounds_v1_0(p1,p2,map));