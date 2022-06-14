function [xx,yy]= return_object_voronoi_points_mod5(L,threshold,otherObjs)
% returnObjectVoronoiPoints takes the nonzero region found in L as returns
% sparse points from boundary and some times cm for - these points can be used
% together with other sets of object-Voronoi-points to make voronoi regions
% (using returnLabeledVoronoi)
%
% The points need to have a certain distance apart (atleast 10 pixels) so that the voronoi
% regions does not become to crowded (we need space for 1-2 pixel wide seperation lines between them)
%% Version history
%  mod3: there were issues with voronoi points not being inside the object;
%  this happens because the points aren't guaranteed integers (CM and
%  extrema) and rounding may lead to a value that is out of the cluster.
%  To address this, will round the resulting points, run a check to see
%  all produced points are inside the cluster; if not, will find nearby
%  one inside the cluster
%  3/27/18: changed interObjTresh to 5 (after some segmentation issues)
%  10/25/19: disconvered an issue with adding in extrema and controids ->
%  the use of '==1' meant often, these points were not added since objects 
%  in L could be any integer. The '==1' was removed so that any integer
%  causes an evaluation of true 
%  mod4: adding a layer of off-boundary voronoi points - hopefully it will
%  make sure object boundaries are better preserved and changed how many
%  boundary points are kept
%  mod5: changing how points are eliminated using distance threshold.
%  Getting rid of the set to (0,0) thing, and increasing the number of
%  boundary points included in the pool

B = bwboundaries(L);
Lslim = bwmorph(L,'thin',threshold+1);
B2 = bwboundaries(Lslim);
% stack boundary points in one n x 2 matrix
nB1 = 0;
nB2 = 0;
for k=1:length(B)
    nB1 = nB1+size(B{k},1);
end
for k=1:length(B2)
    nB2 = nB2+size(B2{k},1);
end
BO = zeros(nB1+nB2,2);
nBk = 0;
for k=1:length(B2)
    BO(nBk+(1:length(B2{k})),:) = B2{k};
    nBk = nBk+length(B2{k});
end
for k=1:length(B)
    BO(nBk+(1:length(B{k})),:) = B{k};
    nBk = nBk+length(B{k});
end

% remove many of boundary points (they are far too close and too many!!!)
cc=0;
freq = round(min(length(BO)*0.5,threshold));
freq = max(1,freq);
temp = zeros(1+floor(round((length(BO)-1)/freq,4)),2);
for k=1:length(BO)
    if mod(k-1,freq)==0
        cc=cc+1;
        temp(cc,:)=BO(k,:);
    end
end

BO = zeros(length(temp),2);

if ~isempty(temp)
    % interchange columns (because bwboudaries returns x and y opposite usual)
    BO(:,1) = temp(:,2);
    BO(:,2) = temp(:,1);
end

% get center of mass and extrema
r = regionprops(logical(L),'Centroid','Extrema');

CM=zeros(length(r),2);
nex = 0;

for k=1:length(r)
    cm = round(r(k).Centroid);
    if(~isempty(cm))
        if(L(cm(2),cm(1)))
            CM(k,:) = cm;
        else
            % try median y
            indlist = find(bwlabel(L)==k);
            [cmy,cmx] = ind2sub(size(L),indlist(1+round(0.5*(length(indlist)-1))));
            cmx2 = round(0.5*(find(L(cmy,1:cmx)==0,1,'last')+find(L(cmy,cmx+1:end)==0,1,'first')+cmx));
            CM(k,:) = [cmx2 cmy];
        end
    end
    % add up number of extrema
    nex = nex + length(r(k).Extrema);
end

EX=zeros(nex,2);
exc = 0;
for k=1:length(r)
    ex = round(r(k).Extrema);
    for kk = 1:size(ex,1)
        if(L(ex(kk,2),ex(kk,1)))
            exc = exc + 1;
            EX(exc,:) = ex(kk,:);
        else
            % if the rounded extrema is not inside the cluster, it should
            % be only 1 pixel away; we will conduct a search
            cut = L(ex(kk,2)+(-1:1),ex(kk,1)+(-1:1));
            [ay,ax] = find(cut);
            if(~isempty(ax))
                dx = ax-2;
                dy = ay-2;
                distsqr = dx.^2+dy.^2;
                [~,imin] = min(distsqr);
                exc = exc + 1;
                EX(exc,:) = ex(kk,:)+[dx(imin) dy(imin)];
            else
                % will skip this extrema
            end
        end
    end
end
EX = EX(1:exc,:);


% add center of masses and extrema otherwise we might loose very small mets or points on the edges of tendrils!!
BO = [CM;EX;BO];
nBO = size(BO,1);

% remove points too close to other objects
num_remaining = nBO;
interObjTresh = min(5,threshold); % 5 is high enough to prevent segmentation errors
was_removed = false(1,nBO);
for k=1:length(otherObjs)
   pointsX = otherObjs(k).voronoiPointsX;
   pointsY = otherObjs(k).voronoiPointsY;
   for n=1:length(pointsX)
       for ii=1:nBO
          dist = sqrt((BO(ii,1)-pointsX(n))^2 + (BO(ii,2)-pointsY(n))^2);
          if (dist < interObjTresh && num_remaining > 1)
              % remove points unless only one point left
              was_removed(ii) = true;
              num_remaining = num_remaining - 1;
          end
       end
   end
end

% if points are too close replace one of them with (0,0) - then later remove
% note earlier points are removed preferentially
for ii=1:nBO
    if(~was_removed(ii))
        for jj=ii+1:nBO
            if(~was_removed(jj))
                dist = sqrt((BO(ii,1)-BO(jj,1))^2 +(BO(ii,2)-BO(jj,2))^2);
                if dist<threshold
                    was_removed(jj) = true;
                end
            end
        end
    end
end

% remove zeros!
xx=BO(~was_removed,1);

yy=BO(~was_removed,2);

end

