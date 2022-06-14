%% findAnchorPoints_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/20/19
%  Project: Fish Image Analysis
%  find anchor points linking endpoints along path

function [anchors] = findAnchorPoints_v1_0(p1,p2,map)

if(isInBounds_v1_0(p1,p2,map))
    anchors = zeros(0,2);
elseif(map(p1(1),p1(2)) && map(p2(1),p2(2)))
    anchor_center = findMidpointAnchor_v1_0(p1,p2,map);
    if(~isempty(anchor_center))
        if(isInBounds_v1_0(p1,anchor_center,map))
            anchor1 = zeros(0,2);
        else
            anchor1 = findAnchorPoints_v1_0(p1,anchor_center,map);
        end
        if(isInBounds_v1_0(anchor_center,p2,map))
            anchor2 = zeros(0,2);
        else
            anchor2 = findAnchorPoints_v1_0(anchor_center,p2,map);
        end
        
        anchors = [anchor1; anchor_center; anchor2];
    else
        warning('No anchor points found between (%i,%i) and (%i,%i). There is a gap in the path.',p1(1),p1(2),p2(1),p2(2));
        anchors = zeros(0,2);
    end
else
    warning('At least one endpoint does not lie on path: (%i,%i)=%i, (%i,%i)=%i.\n Arbitrary anchor points of (1,1) and (1,%i) given',...
        p1(1),p1(2),map(p1(1),p1(2)),p2(1),p2(2),map(p2(1),p2(2)),size(map,2));
    anchors = [1 1; 1 size(map,2)];
end 
end

