function [xx,yy]= return_boundary_voronoi_points(L,threshold)
% returnObjectVoronoiPoints takes the nonzero region found in L as returns
% sparse points from boundary and some times cm for - these points can be used
% together with other sets of object-Voronoi-points to make voronoi regions
% (using returnLabeledVoronoi)
%
% The points need to have a certain distance apart (atleast 10 pixels) so that the voronoi
% regions does not become to crowded (we need space for 1-2 pixel wide seperation lines between them)
%% Version history
%  1.0: modified from return_object_voronoi_points_mod4

B = bwboundaries(L);
% stack boundary points in one n x 2 matrix
nB1 = 0;
for k=1:length(B)
    nB1 = nB1+size(B{k},1);
end
BO = zeros(nB1,2);
nBk = 0;
for k=1:length(B)
    BO(nBk+(1:length(B{k})),:) = B{k};
    nBk = nBk+length(B{k});
end

% remove many of boundary points (they are far too close and too many!!!)
cc=0;
freq = round(threshold*0.5);
freq = max(1,freq);
temp = zeros(1+floor(round((length(BO)-1)/freq,4)),2);
for k=1:length(BO)
    if mod(k-1,freq)==0
        cc=cc+1;
        temp(cc,:)=BO(k,:);
    end
end

BO = zeros(length(temp),2);

if isempty(temp)==0
    % interchange columns (because bwboudaries returns x and y opposite usual)
    BO(:,1) = temp(:,2);
    BO(:,2) = temp(:,1);
end

nBO = size(BO,1);

% if points are too close replace one of them with (0,0) - then later remove
% note earlier points are removed preferentially
was_removed = false(1,nBO);
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
xx = BO(~was_removed,1);
yy = BO(~was_removed,2);

end

