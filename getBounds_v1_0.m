%% getBounds_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/19/19
%  Project: Fish Image Analysis
%  see which pixels in line of pixels in between two points lies on path
function [bounds,line] = getBounds_v1_0(p1,p2,map)

% fprintf('(%i,%i) -> (%i,%i)\n',p1(1),p1(2),p2(1),p2(2));
line = getPixelLine_v1_0(p1,p2);
nrow = size(map,1);
ii = line(:,1)+(line(:,2)-1)*nrow;
bounds = map(ii);