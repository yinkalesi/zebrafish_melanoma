function bw  = smoothBW_mod(bw,delta)
% smooths the outline of connected regions of 1's in images with only 0's
% and 1's

B = bwboundaries(bw);
L = bwlabel(bw);

temp = false(size(L));

for k=1:length(B)
boundary = B{k};
% smooth_boundary gives an error if fewer boundary points than delta - the
% fix is that delta is set to number of boundary points for small sets - 
% boundary segments smaller than delta will be smoothed down to one point
used_delta = min(length(boundary(:,1)),delta);
[boundary(:,2), boundary(:,1)] = smooth_boundary(boundary(:,1), boundary(:,2),used_delta);

    for i=1:length(boundary)
        temp(round(boundary(i,2)),round(boundary(i,1)))=1;
    end

end

temp=imfill(temp,'holes');

bw=temp;

