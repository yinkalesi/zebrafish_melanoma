function [bw_edge,rect_tail] = getTailEdgeImage_v1_0(bw,bf,tail_coor)
%% getTailEdgeImage_v1_0
%  Version 1.0
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 10/21/17
%  Project: Tumor Growth, Logarithmic Continuum Form

% crop image to tail region
% right image
bws = bw;
tail_wid = mean(sum(bws(:,tail_coor(1)-[150 50])));
rect_tail = round([tail_coor(1)-50 tail_coor(2)-0.5*tail_wid 60 tail_wid]);
tail = imcrop(bf,rect_tail);
% brighten image (looking for a darkline so multiplying the intensities
% amplifies the brightness of the bright spots compared to the dark
amp = min(3,65535/mean2(tail)*0.8);
tail = tail*amp; % note tail is uint16 and has a max value of 65535; 
% this operation essentially washes out any part of tail image higher than average intensity

% create cropped mask of body for tail region
bwc = imcrop(bws,rect_tail);
% assuming the line we want is darker than the average pixel, we can create
% a mask to screen for dark pixels
tail_int = double(tail)/double(max(tail(:)));
dark_thresh = 0.9;
bw_dark = 1-imbinarize(tail_int,dark_thresh);
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

% figure(2);
% subplot(1,2,1);
% imshow(bw_edge);
% subplot(1,2,2);
% imshow(double(tail).*double(bwc),[])
% pause(0.1);
end