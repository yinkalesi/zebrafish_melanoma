function ch = getConvHull(bw)
%% getConvHull
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 11/10/17
% get convex hull image from bw image

[ys,xs] = find(bw==1);
k = convhull(xs,ys);
ch = poly2mask(xs(k),ys(k),size(bw,1),size(bw,2));
ch = bwmorph(ch,'thicken',1);