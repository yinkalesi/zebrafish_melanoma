function bw = return_eye_outline_mod4(bf_in,rect)
% takes rect and red channel of rgb bf image and returns eye outline
% mod2: added blurring
% blur bf image to get more clear-cut results
% mod3: unsure circularity using max height as a guide for diameter
% mod4: refinements

xblur = [1 1 2 1 1];
yblur = [1 1 2 1 1]';
blur = yblur*xblur;
blur = blur/sum(blur(:));
bf = imfilter(double(bf_in)/65535,blur,'replicate');
crop = imcrop(bf,rect);
% don't want crop to include any really bright regions becuause
% this could mess with thresholding
crop_list = crop(crop<5e4/65535);

bw = imbinarize(bf,graythresh(crop_list));
bw = bw.*-1 +1;

bw(:,1:round(rect(1)/2)) = 0;
[row,col] = size(bw);
% have to be careful with crop in case rect doesn't lie where it should on
% the eye
% expected_eye_start = rect(1)-rect(3)/10;
first_eye_locs = find(bw(:),10)/row;
first_eye_col = mean(first_eye_locs);
actual_eye_start = floor(first_eye_col);
eye_center = round(mean(first_eye_locs-floor(first_eye_locs))*row);
erasure_start = max(actual_eye_start+100,rect(1)+round(rect(3)*.64));
% remove large circle around eye
xdv = (1:col)-0.5*(actual_eye_start+erasure_start);
ydv = (1:row)-eye_center;
[xds,yds] = meshgrid(xdv.^2,ydv.^2);
dist_from_eye_sqr = xds+yds;
bw_circ100 = dist_from_eye_sqr <= 0.25*(erasure_start-actual_eye_start)^2;
bw(bw_circ100==0) = 0;

bw = bwmorph(bw,'close',2); % get rid of isolated false positive regions
bw = bwmorph(bw,'open',2);
bw = bwmorph(bw,'dilate',2);

bw = imfill(bw,'holes');

bw = bwareaopen(bw,300);

bw = bwmorph(bw,'erode',2);
bw = bwareaopen(bw,300);

bw = bwmorph(bw,'dilate',2);

% there should only be one large object in image - just incase,
% will try and remove any moderately large objects
area = length(find(bw>0));
bw = bwareaopen(bw,round(area*0.25));
bw = smoothBW_mod(bw,10);

% 1) height of eye area - should fit h^2 = D^2-4*(c-x)^2
% where D is the max height and c = x value at max height
xstart = floor(find(bw,1,'first')/size(bw,1))+1;
xend = floor(find(bw,1,'last')/size(bw,1))+1;
h = zeros(1,xend-xstart+1);
cyh = zeros(1,xend-xstart+1);
for i = 1:xend-xstart+1
    line = bw(:,i+xstart-1);
    top = find(line,1,'first');
    bot = find(line,1,'last');
    if(~isempty(top))
        h(i) = bot-top+1;
        cyh(i) = 0.5*(top+bot);
    end
end
% 2) find center and diameter of eye
if(~isempty(h))
    ihmax = round(2/3*length(h));
    [D,iD] = max(h(1:ihmax));
    cx = 0.5*D+xstart-1;
    cy = cyh(iD);
    % 3) remove anything in bw outside of circle on the right side
    xdv = (1:col)-cx;
    ydv = (1:row)-cy;
    [xds,yds] = meshgrid(xdv.^2,ydv.^2);
    dist_from_eye_sqr = xds+yds;
    bw_circ = dist_from_eye_sqr <= 0.25*D^2;
    bw_circ(:,1:floor(cx)) = 1;
    bw(bw_circ==0) = 0;
end
end

