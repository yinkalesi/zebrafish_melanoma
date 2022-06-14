function bw  = smoothBW(bw,delta)
% smooths the outline of connected regions of 1's in images with only 0's
% and 1's

B = bwboundaries(bw);
L = bwlabel(bw);

temp = zeros(size(L));

for k=1:length(B)
boundary = B{k};

[boundary(:,2), boundary(:,1)] = smooth_boundary(boundary(:,1), boundary(:,2),delta);

    for i=1:length(boundary)
        temp(round(boundary(i,2)),round(boundary(i,1)))=1;
    end

end

temp=imfill(temp,'holes');

bw=temp;

