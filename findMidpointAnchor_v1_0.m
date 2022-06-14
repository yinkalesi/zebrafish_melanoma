%% findMidpointAnchor_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/19/19
%  Project: Fish Image Analysis
%  find the closest point on tangent line across midpoint on the path
function [anchor] = findMidpointAnchor_v1_0(p1,p2,map)

msize = size(map);
[~,line] = getBounds_v1_0(p1,p2,map);
midpoint = line(round(length(line)*0.5),:);
tanslope = -(p2(1)-p1(1))/(p2(2)-p1(2));

% find two edge points along tangent
y_x1 = round(midpoint(2)+tanslope*(1-midpoint(1)));
if y_x1 > msize(2)
    tp1 = round([(msize(2)-midpoint(2))/tanslope+midpoint(1) msize(2)]);
elseif y_x1 < 1
    tp1 = round([(1-midpoint(2))/tanslope+midpoint(1) 1]);
else
    tp1 = round([1 y_x1]);
end

y_xe = round(midpoint(2)+tanslope*(msize(1)-midpoint(1)));
if y_xe > msize(2)
    tp2 = round([(msize(2)-midpoint(2))/tanslope+midpoint(1) msize(2)]);
elseif y_xe < 1
    tp2 = round([(1-midpoint(2))/tanslope+midpoint(1) 1]);
else
    tp2 = round([msize(1) y_xe]);
end
    
% find an intersect
[inter1,dist1] = findIntersect_v1_0(midpoint,tp1,map);
[inter2,dist2] = findIntersect_v1_0(midpoint,tp2,map);

if(isempty(inter1) && ~isempty(inter2))
    anchor = inter2;
elseif(isempty(inter2) && ~isempty(inter1))
    anchor = inter1;
elseif(~isempty(inter1) && ~isempty(inter2))
    if(dist1<=dist2)
        anchor = inter1;
    else
        anchor = inter2;
    end
else
    % warning('No viable midpoint between (%i,%i) and (%i,%i)',p1(1),p1(2),p2(1),p2(2));
    anchor = [];
end
% disp(anchor);
    


