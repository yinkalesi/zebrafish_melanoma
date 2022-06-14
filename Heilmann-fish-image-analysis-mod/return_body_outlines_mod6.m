function [body1,body2,body1_end,body2_end,edge1,edge2,edge_box1,edge_box2,left_shift] = return_body_outlines_mod6(bf1,R_body1,bf2,R_body2,min_con_reg)
% Versions:
% mod: bwareaopen threshold changed from 20000 to 0.5*min_con_reg
% mod2: change out fins are removed and tail point is set
% mod3: return tail edge image
% mod4: better tail identification
% mod5: avoids large alignment adjustments, 
%       change test_depth/average intensity testing approach
% mod6: trying to fix issue in one fish

% erode body (find 'core' of tail - get rid of tail fin)
S_body1 = bwmorph(R_body1,'erode',5);
% S_body1 = bwmorph(S_body1,'thin',20);
S_body1 = bwareaopen(S_body1,0.5*min_con_reg);
S_body2 = bwmorph(R_body2,'erode',5);
% S_body2 = bwmorph(S_body2,'thin',20);
S_body2 = bwareaopen(S_body2,0.5*min_con_reg);

% find nose and tail coor of both left and right side images
 
[y,x] = find(S_body1,10,'first');
nose_coor1 = [round(mean(x)),round(mean(y))];
[y,x] = find(S_body1,100,'last');
tail_coor1 = [round(mean(x)),round(mean(y))];

[y,x] = find(S_body2,10,'first');
nose_coor2 = [round(mean(x)),round(mean(y))];
[y,x] = find(S_body2,100,'last');
tail_coor2 = [round(mean(x)),round(mean(y))];

% look for better tail_coor (sometimes the rough outline doesn't cover the
% tail)
% detect if long section fin was not removed well enough by looking at
% intensity at end of body
test_depth = 200;
S_body1_long = S_body1;
S_body1_long(tail_coor1(2)+(-20:20),tail_coor1(1)-5:min(tail_coor1(1)+test_depth,end)) = 1;
S_body2_long = S_body2;
S_body2_long(tail_coor2(2)+(-20:20),tail_coor2(1)-5:min(tail_coor2(1)+test_depth,end)) = 1;
bf1s = bf1(:,tail_coor1(1)-test_depth+1:min(tail_coor1(1)+test_depth,end));
bf1s(S_body1_long(:,tail_coor1(1)-test_depth+1:min(tail_coor1(1)+test_depth,end))==0) = 0;
sbf1s = sum(bf1s,1);
ssb1 = sum(S_body1_long(:,tail_coor1(1)-test_depth+1:min(tail_coor1(1)+test_depth,end)),1);
abf1s = sbf1s./(ssb1+0.001);
bf2s = bf2(:,tail_coor2(1)-test_depth+1:min(tail_coor2(1)+test_depth,end));
bf2s(S_body2_long(:,tail_coor2(1)-test_depth+1:min(tail_coor2(1)+test_depth,end))==0) = 0;
sbf2s = sum(bf2s,1);
ssb2 = sum(S_body2_long(:,tail_coor2(1)-test_depth+1:min(tail_coor2(1)+test_depth,end)),1);
abf2s = sbf2s./(ssb2+0.001);
% propose new end point and test for acceptance
prop_sep1 = getSpatialThreshold_v1_2(abf1s);
left_avg1 = mean(abf1s(1:prop_sep1-1));
right_avg1 = mean(abf1s(prop_sep1:end));
prop_sep2 = getSpatialThreshold_v1_2(abf2s);
left_avg2 = mean(abf2s(1:prop_sep2-1));
right_avg2 = mean(abf2s(prop_sep2:end));
orig_tailx1 = tail_coor1(1);
orig_tailx2 = tail_coor2(1);
if(right_avg1>left_avg1) % && left_avg1/right_avg1 < 0.6667)
    tail_coor1(1) = tail_coor1(1)-test_depth+prop_sep1;
end
if(right_avg2>left_avg2) % && left_avg2/right_avg2 < 0.6667)
    tail_coor2(1) = tail_coor2(1)-test_depth+prop_sep2;
end
% use cropping to get better identification of features
rough_tail1 = imcrop(bf1(:,:,1),[tail_coor1(1)-100 tail_coor1(2)-50 300 100]);
rough_tail2 = imcrop(bf2(:,:,1),[tail_coor2(1)-100 tail_coor2(2)-50 300 100]);
% tail should be lighter than body - the end point of body should be when
% the brightest sharply increases (transition from body to tail)
int_rt1 = double(mean(rough_tail1,1))/double(max(rough_tail1(:)));
int_rt2 = double(mean(rough_tail2,1))/double(max(rough_tail2(:)));
% thresh_irt1 = graythresh(int_rt1);
% thresh_irt2 = graythresh(int_rt2);
% better_tailx1 = find(int_rt1<thresh_irt1,1,'last')+tail_coor1(1)-101;
% better_tailx2 = find(int_rt2<thresh_irt2,1,'last')+tail_coor2(1)-101;
better_tailx1 = getSpatialThreshold_v1_2(int_rt1)+tail_coor1(1)-101;
better_tailx2 = getSpatialThreshold_v1_2(int_rt2)+tail_coor2(1)-101;
% old_tailx1 = tail_coor1(1);
% old_tailx2 = tail_coor2(1);
tail_coor1(1) = better_tailx1;
tail_coor2(1) = better_tailx2;
% calculate fish lengths
fish_length1 = tail_coor1(1)-nose_coor1(1);
fish_length2 = tail_coor2(1)-nose_coor2(1);
% rectangle for looking at only end of tail (this is where threshold for finding boundary is hardest to determine)
rect_tail1 = [nose_coor1(1) + fish_length1.*0.95  tail_coor1(2)-100 200 200];
rect_tail2 = [nose_coor2(1) + fish_length2.*0.95  tail_coor2(2)-100 200 200];

% cut out end of tail in red channel of rgb image - use this image
% to find theshold
tail_im1 = imcrop(bf1(:,:,1),rect_tail1);
tail_im2 = imcrop(bf2(:,:,1),rect_tail2);
tail_rb1 = imcrop(R_body1,rect_tail1);
tail_rb2 = imcrop(R_body2,rect_tail2);

% find theshold using graythres (uses Otsu's method) - the brightest pixels
% will not be supplied to the thresholding algorithm to hopefully get
% better separation between the fish body and the tail fin
mean_bright1 = mean(tail_im1(tail_rb1==0));
mean_bright2 = mean(tail_im2(tail_rb2==0));
max_bf1 = 65535;
max_bf2 = 65535;
tail_list1 = double(tail_im1(tail_im1<mean_bright1))/mean_bright1;
tail_list2 = double(tail_im2(tail_im2<mean_bright2))/mean_bright2;
T1 = graythresh(tail_list1)*mean_bright1/max_bf1;
T2 = graythresh(tail_list2)*mean_bright2/max_bf2;

% Make bw images using fish specific thresholds
bw1 = 1-imbinarize(bf1(:,:,1),T1);
bw2 = 1-imbinarize(bf2(:,:,1),T2);

% need to make sure S_body extends past tail
S_body1(tail_coor1(2)+(-5:5),orig_tailx1:better_tailx1)=1;
S_body2(tail_coor2(2)+(-5:5),orig_tailx2:better_tailx2)=1;
% dilate rough body outline to use as mask for new body outline (remove dirt and shadows around fish)
S_body1 = bwmorph(S_body1,'dilate',38);
S_body1 = bwareaopen(S_body1,0.5*min_con_reg);
S_body2 = bwmorph(S_body2,'dilate',38);
S_body2 = bwareaopen(S_body2,0.5*min_con_reg);
% set head region part to be same as R_body (R_body does a decent job of
% identifying the head region)
S_body1(:,1:round(nose_coor1(1)+0.25*fish_length1)) = R_body1(:,1:round(nose_coor1(1)+0.25*fish_length1));
S_body2(:,1:round(nose_coor2(1)+0.25*fish_length2)) = R_body2(:,1:round(nose_coor2(1)+0.25*fish_length2));
S_body1 = bwmorph(S_body1,'dilate',2);
S_body2 = bwmorph(S_body2,'dilate',2);
% use mask
bw1(S_body1==0)=0;
bw2(S_body2==0)=0;
bw1(:,round(nose_coor1(1) + fish_length1.*1.2):end)=0;
bw2(:,round(nose_coor1(1) + fish_length2.*1.2):end)=0;

% fill hole, erode and dilate to get rid of small protrusions
% (smooths boundary)
% tbws are cleaned up more harshly (there may be holes in fish body)

tbw1 = bwmorph(bw1,'erode',10);
tbw2 = bwmorph(bw2,'erode',10);
tbw1 = bwareaopen(tbw1,0.5*min_con_reg);
tbw2 = bwareaopen(tbw2,0.5*min_con_reg);
tbw1 = bwmorph(tbw1,'dilate',11);
tbw2 = bwmorph(tbw2,'dilate',11);
% use tbws to get convex hull that (hopefully) doesn't contain any fins
ch1 = getConvHull(tbw1);
ch2 = getConvHull(tbw2);
ch1(:,[1:round(nose_coor1(1)+0.25*fish_length1) round(tail_coor1(1)-0.25*fish_length1):end])=1;
ch2(:,[1:round(nose_coor2(1)+0.25*fish_length2) round(tail_coor2(1)-0.25*fish_length2):end])=1;
% clean up bws gently
bw1 = bwfill(bw1,'holes');
bw2 = bwfill(bw2,'holes');
bw1 = bwmorph(bw1,'thicken',9);
bw2 = bwmorph(bw2,'thicken',9);
bw1 = bwmorph(bw1,'dilate',1);
bw2 = bwmorph(bw2,'dilate',1);
bw1 = bwmorph(bw1,'thin',9);
bw2 = bwmorph(bw2,'thin',9);
bw1 = bwmorph(bw1,'erode',1);
bw2 = bwmorph(bw2,'erode',1);
bw1 = bwareaopen(bw1,0.5*min_con_reg); % remove regions of 1's smaller than 0.5*min_con_reg
bw2 = bwareaopen(bw2,0.5*min_con_reg); % remove regions of 1's smaller than 0.5*min_con_reg
bw1 = imfill(bw1,'holes');
bw2 = imfill(bw2,'holes');
% use chs to try and get rid of fins
bw1(ch1==0) = 0;
bw2(ch2==0) = 0;

% need to find better approximation for end of tail
% there are vertical lines signaling the beginning of the tail - we will
% try to identify this line

% it is useful to get rid of any thin tail fin parts that have not been
% removed already
% as of mod4, already used brightness method above to improve tail fin
% identification and the tail width method below seems less reliable than
% the brightness method - set tail_cut to a smaller number
tail_cut = 50;
endrow1 = min(tail_coor1(1)+tail_cut,size(bw1,2));
endrow2 = min(tail_coor2(1)+tail_cut,size(bw2,2));
% mean_cut = round(tail_cut*0.5);
tail_wid1 = sum(bw1(:,tail_coor1(1)-tail_cut:endrow1));
tail_wid2 = sum(bw2(:,tail_coor2(1)-tail_cut:endrow2));
% mean_wid1 = mean(tail_wid1(1:mean_cut));
% mean_wid2 = mean(tail_wid2(1:mean_cut));
% otsu_wid1 = graythresh(tail_wid1/max(tail_wid1))*max(tail_wid1);
% otsu_wid2 = graythresh(tail_wid2/max(tail_wid2))*max(tail_wid2);
% thresh_wid1 = max(mean_wid1*0.5,otsu_wid1);
% thresh_wid2 = max(mean_wid2*0.5,otsu_wid2);
% tail_end1 = tail_coor1(1)-tail_cut-1+find(tail_wid1>thresh_wid1,1,'last');
% tail_end2 = tail_coor2(1)-tail_cut-1+find(tail_wid2>thresh_wid2,1,'last');
izero1 = find(tail_wid1==0,1,'first');
if(isempty(izero1))
    izero1 = length(tail_wid1);
end
izero2 = find(tail_wid2==0,1,'first');
if(isempty(izero2))
    izero2 = length(tail_wid2);
end
tail_end1 = getSpatialThreshold_v1_2(tail_wid1(1:izero1))+tail_coor1(1)-tail_cut-1;
tail_end2 = getSpatialThreshold_v1_2(tail_wid2(1:izero2))+tail_coor2(1)-tail_cut-1;
% cut of tail fin area
bw1(:,max(tail_end1,10+tail_coor1(1)):end)=0;
bw2(:,max(tail_end2,10+tail_coor2(1)):end)=0;

% get image with better tail region
[body1_end_ft,tail1,edge1,edge_box1] = getFlatTailEnd(bw1,bf1(:,:,1));
[body2_end_ft,tail2,edge2,edge_box2] = getFlatTailEnd(bw2,bf2(:,:,1));

% align edges to get better alignment for matching both sides of fish
[edge1b,dr,dc] = alignBW_v2_0(edge1,edge2,1,1,10,10);

% if dr or dc is too high
dlim = 20;
if(abs(dr) > dlim)
    warning('set dr to 0 because shift was too large');
    dr = 0;
end
if(abs(dc) > dlim)
    warning('set dc to 0 because shift was too large');
    dc = 0;
end

% use alignment to get final body end
% calculate relative positions of end points inside cropped image
rel_pos1 = body1_end_ft(1)-edge_box1(1)+1;
rel_pos2 = body2_end_ft(1)-edge_box2(1)+1;
% select relative position desired for fish tail end position out of the
% two available from the two sides of the fish
center_pos = min(edge_box1(3),edge_box2(3))*0.66;
if(abs(rel_pos1-center_pos) < abs(rel_pos2-center_pos))
    selected_pos1 = rel_pos1;
else
    selected_pos1 = rel_pos2;
end
% the alignment gives states how many pixels edge1 needs to be shifted to
% match edge2 - this means the coordinates of the crop rectangle can used
% to relate positions on one side of the fish to the other
body1_end(1) = edge_box1(1)+selected_pos1-1-dc;
body2_end(1) = edge_box2(1)+selected_pos1-1;
selected_pos2 = min(edge_box1(4)*0.5,edge_box2(4)*0.5);
body1_end(2) = edge_box1(2)+selected_pos2-1-dr;
body2_end(2) = edge_box2(2)+selected_pos2-1;
% position of edge image relative to body image
left_shift = [dr dc];

% modify body images
body1 = bw1;
body1(:,body1_end(1):end) = 0;
body2 = bw2;
body2(:,body2_end(1):end) = 0;

% make boundary of fish body smoother
body1 = smoothBW_mod(body1,10);
body2 = smoothBW_mod(body2,10);

% % make plots
% figure(2);
% subplot(2,3,1);
% imshow(tail1);
% hold on;
% plot(body1_end(1)-edge_box1(1)+1,body1_end(2)-edge_box1(2)+1,'o');
% plot(body1_end_ft(1)-edge_box1(1)+1,body1_end_ft(2)-edge_box1(2)+1,'*');
% plot(body1_end(1)*[1 1]-edge_box1(1)+1,[1 edge_box1(4)+1]);
% hold off;
% title(sprintf('Tail end: %i',body1_end(1)));
% subplot(2,3,4);
% imshow(tail2);
% hold on;
% plot(body2_end(1)-edge_box2(1)+1,body2_end(2)-edge_box2(2)+1,'o');
% plot(body2_end_ft(1)-edge_box2(1)+1,body2_end_ft(2)-edge_box2(2)+1,'*');
% plot(body2_end(1)*[1 1]-edge_box2(1)+1,[1 edge_box2(4)+1]);
% hold off;
% title(sprintf('Tail end: %i',body2_end(1)));
% subplot(2,3,2);
% imshow(edge1);
% hold on;
% plot(body1_end(1)-edge_box1(1)+1,body1_end(2)-edge_box1(2)+1,'o');
% plot(body1_end_ft(1)-edge_box1(1)+1,body1_end_ft(2)-edge_box1(2)+1,'*');
% plot(body1_end(1)*[1 1]-edge_box1(1)+1,[1 edge_box1(4)+1]);
% hold off;
% subplot(2,3,5);
% imshow(edge2);
% hold on;
% plot(body2_end(1)-edge_box2(1)+1,body2_end(2)-edge_box2(2)+1,'o');
% plot(body2_end_ft(1)-edge_box2(1)+1,body2_end_ft(2)-edge_box2(2)+1,'*');
% plot(body2_end(1)*[1 1]-edge_box2(1)+1,[1 edge_box2(4)+1]);
% hold off;
% subplot(2,3,3)
% imshowpair(getShiftedImg_v1_0(tail1,dr,dc),tail2);
% title(sprintf('dr: %i,dc: %i',dr,dc));
% subplot(2,3,6)
% imshowpair(edge1b,edge2);
% hold on;
% plot(body2_end(1)*[1 1]-edge_box2(1)+1,[1 edge_box2(4)+1]);
% hold off;
% pause(0.1);
% % if(abs(dr)>50 || abs(dc)>50)
% %     disp('Stop');
% % end
end

function [tail_end,tail_int,bw_edge,rect_tail] = getFlatTailEnd(bw,bf)
% crop image to tail region
% right image
bws = smoothBW_mod(bw,5);
[tail_coor, ~] = getBodyEndPoint_v1_0(bws,50,150);
tail_wid = mean(sum(bws(:,tail_coor(1)-[150 50])));
% [y,x] = find(bws,1000,'last');
% tail_coor = [round(mean(x)),round(mean(y))];
rect_tail = round([tail_coor(1)-65 tail_coor(2)-0.5*tail_wid 65 tail_wid]);
tail = imcrop(bf,rect_tail);
% brighten image (looking for a darkline so multiplying the intensities
% amplifies the brightness of the bright spots compared to the dark
amp = min(3,65535/mean2(tail)*0.8);
tail = tail*amp; % note tail is uint16 and has a max value of 65535; 
% this operation essentially washes out any part of tail image higher than average intensity

% create cropped mask of body for tail region
bws = bwmorph(bws,'erode',5);
bwc = imcrop(bws,rect_tail);
% assuming the line we want is darker than the average pixel, we can create
% a mask to screen for dark pixels
tail_int = double(tail)/double(max(tail(:)));
dark_thresh = 0.9;
bw_dark = imbinarize(tail_int,dark_thresh)*-1+1;
bw_dark = bwmorph(bw_dark,'thicken',1);
bw_dark(bwc==0)=0;
% get gradient image 
hx = fspecial('sobel')';
hx = hx/sum(abs(hx(:)))*2; % maximum grad value magnitude is 1
xweights = [1 4 1];
xweights = xweights/sum(xweights);
yweights = [1 2 4 2 1]';
yweights = yweights/sum(yweights);
blur = yweights*xweights;
grad = imfilter(tail_int,hx,'replicate');
grad = imfilter(grad,blur,'replicate'); % a vertical blur effect
grad(bw_dark==0)=0;
% split positive and negative
pos_grad = max(0,grad);
neg_grad = max(0,-grad);
% want to remove noise do to random changes in fish color
% will use two thresholding approaches. Otsu's method (graythresh) usually
% works but fails in at least one known case: when a highly pigmented tumor
% appears close to the tail, resulting in a dark to light transition that
% creates high values in neg_grad and making the calculated threshold high;
% an alternative method using the mean and std of the data is added to
% avoid this
p_thresh1 = graythresh(nonzeros(pos_grad(:)));
[p_avg,p_std] = getBackgroundMean_v1_0(pos_grad);
p_thresh2 = p_avg+1.5*p_std;
n_thresh1 = graythresh(nonzeros(neg_grad(:)));
[n_avg,n_std] = getBackgroundMean_v1_0(neg_grad);
n_thresh2 = n_avg+1.5*n_std;
pos_thresh = min(p_thresh1,p_thresh2);
neg_thresh = min(n_thresh1,n_thresh2);
pos_grad(pos_grad<pos_thresh)=0;
neg_grad(neg_grad<neg_thresh)=0;
% bwpg = imbinarize(pos_grad,graythresh((pos_grad(:))));
bw_fill = logical(fillHorizontalEdges_v1_0(-pos_grad,neg_grad));
bw_edge = bwmorph(bw_fill,'spur',Inf);
% bridge 1 pixel vertical gaps:
bw_edge = imfilter(double(bw_edge),[1 1 1]'/3);
bw_edge(bw_edge<0.5) = 0;
bw_edge = logical(bw_edge);
% bw_edge = bwmorph(bw_edge,'thin',2);
% bw_edge = bwmorph(bw_edge,'spur',4);
% bwpg2 = imfilter(double(bwedge),[1 0 1]);
% bwedge(bwpg2==2)=0;
bw_edge = bwmorph(bw_edge,'clean');

% find tallest object - this is expected to be the tail end
locs = regionprops(bw_edge,'Centroid','Extrema');
lens = table2array(regionprops('table',bw_edge,'MajorAxisLength'));
orien1 = table2array(regionprops('table',bw_edge,'Orientation'));
% only want objects that are reasonable vertial (60-90 deg)
sellist  = find(90-abs(orien1)<45);
[~,itall] = max(lens(sellist));
itall = sellist(itall);
% make tallest object x position the tail of fish
if(~isempty(itall))
    centroid = locs(itall).Centroid;
    [~,iext] = max(centroid(:,1));
    xtall = round(centroid(iext,1));
    tail_end = [rect_tail(1)+xtall-1,rect_tail(2)+round(centroid(iext,2))-1];
elseif(~isempty(locs))
    [~,itall] = max(lens);
    centroid = locs(itall).Centroid;
    [~,iext] = max(centroid(:,1));
    xtall = round(centroid(iext,1));
    tail_end = [rect_tail(1)+xtall-1,rect_tail(2)+round(centroid(iext,2))-1];
else
    % no change if tail identification doesn't work
    warning(['Failed tail identification in ' mfilename]);
    tail_end = [tail_coor(1)-30,tail_coor(2)];
end

% % modify bw images
% bw_out = bw;
% bw_out(:,new_tail(1):end) = 0;
% 
% % make boundary of fish body smoother
% bw_out = smoothBW_mod(bw_out,10);

% figure(2);
% subplot(1,2,1);
% imshow(bw_edge);
% hold on;
% if(~isempty(itall))
%     plot(centroid(iext,1),centroid(iext,2),'*');
% end
% plot((tail_end(1)-rect_tail(1)+1)*[1 1],[1 rect_tail(4)+1]);
% hold off;
% subplot(1,2,2);
% imshow(double(tail).*double(bwc),[])
% hold on;
% if(~isempty(itall))
%     plot(centroid(iext,1),centroid(iext,2),'*');
% end
% plot((tail_end(1)-rect_tail(1)+1)*[1 1],[1 rect_tail(4)+1]);
% hold off;
% pause(0.1);
end