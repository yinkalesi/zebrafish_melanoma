%% return_side_mod13
%  Version 13
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 2/10/20
%  Project: Tumor Growth, Logarithmic Continuum Form
% flourescentImage is type uint16
%% Version History
% mod2: putting back varying threshold detection
% mod4: add arguments for threshold minimums
% mod6: really trying to do tumor detection better - the key issue is that
% the thresholding doesn't always work: different tumors have different
% brightness and there is sometimes a halo effect around bright tumors - a
% high threshold misses tumors but a low tumor puts some of the halo region
% into the tumor area. It seems like edge detection is a better paradigm to
% detect tumors. Edge detection was used in previous versions but it wasn't
% the primary method.
% mod7: issue arose with detecting a signal when the background is all
% zeros (since zeros are removed to do the analysis, all the analyzed
% pixels turn out to be bright signal pixels and attempts to separate out
% background pixels result in removing good data)
% mod8: more parameter inputs so sensitivity can be adjusted to process rfp
% images
% mod9: using a low threshold to screen
% mod10: change rescaling of fluorescentImage, using upper threshold
% mod11: It's 2019. Let's change it up. Trying to sovle - 1) False
% negatives where dimmer clusters aren't picked up when bright clusters are
% around; 2) False positives when everything is dim and some random region
% is selected as a tumor
% mod12: problem of losing track of shrinking tumors early, before they
% truly disappear

function [bw, candidates] = return_side_mod13(fluorescentImage,body,eye,specs,brinum,edgefac,ref_int,past_tumors)

if(~exist('past_tumors','var'))
    past_tumors = false(size(fluorescentImage));
end

eye = logical(eye);

% make into type double with value between 0 and 1
fp = min(1,double(fluorescentImage)./double(ref_int));
fp(body==0)=0; % want to concentrate on tumor signal only

% brightest cutoff set by user
thresh_highest = specs.bw_threshold_high*65535/double(ref_int);

% cummulative distribution of intensities
fp_list = fp(body~=0);
nbin = 1000;
[cdf,bin] = histcounts(fp_list,nbin,'Normalization','cdf');
pdf = (cdf-[0 cdf(1:end-1)])*nbin;
xpdf = 0.5*(bin(1:end-1)+bin(2:end));
smooth_pdf = smoothingCurve_v1_0(pdf,10);
smoother_pdf = smoothingCurve_v1_0(pdf,100);
adj_pdf = max(0,smooth_pdf-smoother_pdf);
% dpdf = ([smooth_pdf(2:end) 0]-[0 smooth_pdf(1:end-1)])*nbin;

% find brightest region (brightest 1000 pixels)
body_area = sum(body(:));
per_bnum = brinum/body_area;
i_bnum = find(cdf>=1-per_bnum,1);
thresh_bnum = bin(i_bnum);
% find the 95th percentile of nonzero pixels
nzfp = nonzeros(fp_list);
num_nzfp = length(nzfp);
per5nz = 0.05*num_nzfp/body_area;
i5nz = find(cdf>=1-per5nz,1);
thresh5nz = bin(i5nz);

% use Otsu's method to identify candidates - the issue with Otsu's method
% is that the cutoff can be too high or too low depending on various
% factors but it usually works well
thresh_otsu = graythresh(nzfp);

% find brightest peak
peak_th = brinum/body_area/.05;
peak_cands = find(adj_pdf>=peak_th&xpdf>=thresh_otsu);
if(isempty(peak_cands))
    thresh_brightest = max(thresh_bnum,thresh5nz);
else
    peak_cand_gaps = peak_cands(2:end)-peak_cands(1:end-1);
    last_peak_start = find(peak_cand_gaps>1,1,'last');
    if(isempty(last_peak_start))
        last_peak_start = 1;
    end
    thresh_brightest = xpdf(peak_cands(last_peak_start));
end

bw_brightest = imbinarize(fp,thresh_brightest);
bw_brightest = bwmorph(bw_brightest,'clean');



% get background intensity
[avg_fp,std_fp] = getBackgroundMean(fp_list);
[avg_fp2,std_fp2] = getBackgroundMedian(fp_list);
% want to include more pixels as number of nonzero pixels gets smaller -
% this means the variance between nz pixels is higher
sepfac = 3*max(0,1-3*std_fp2);
thresh_signal = min(avg_fp+sepfac*std_fp,avg_fp2+sepfac*std_fp2);
% signal significantly above background (3 std) should be investigated
bw_signal = imbinarize(fp,thresh_signal);
bw_signal = bwareaopen(bw_signal,3);


% sometimes Otsu's method threshold is much too low (usually in rfp images
% with few bright spots) - this results in >50 of the body being selected
if(thresh_otsu>min(avg_fp2+0.5*std_fp,avg_fp+0.5*std_fp2))
    bw_otsu = imbinarize(fp,thresh_otsu);
    bw_otsu = bwareaopen(bw_otsu,3);
else
    bw_otsu = false(size(fp));
end

% absolute minimum signal that can be accepted as tumor (currently only
% applies for tumors detected solely by edge detection
min_int_thresh = 0.0039; % hard coding this right now - will try and optimize later
bw_min = imbinarize(fp,min_int_thresh*65535/double(ref_int));
bw_min = bwmorph(bw_min,'dilate',1);
bw_min = bwmorph(bw_min,'thin',1);

% intensity gradient
sobely = 0.25*fspecial('sobel');
sobelx = sobely';
gx = imfilter(fp,sobelx,'replicate');
thinbody = bwmorph(body,'erode',2);
gx(thinbody==0)=0;
gy = imfilter(fp,sobely,'replicate');
gy(thinbody==0)=0;
% separate positive and negative, and use blur filter
xweights = [1 4 6 4 1];
xweights = xweights/sum(xweights);
yweights = [1 4 6 4 1]';
yweights = yweights/sum(yweights);
blurx = yweights*xweights;
blury = blurx';

gxp = max(0,imfilter(gx,blurx,'replicate'));
gxn = max(0,imfilter(-gx,blurx,'replicate'));
gyp = max(0,imfilter(gy,blury,'replicate'));
gyn = max(0,imfilter(-gy,blury,'replicate'));

% clean up noise in edge images; it is reasonable to assume the edge pixels
% are rare compared to background pixels so that the edge pixels don't
% affect the mean and std of the image very much. Therefore, we can assume
% pixels within 2 std of the mean are background pixels
% Alternatively, will use the median to get an 'average' (In hopes of
% getting better results for when there is not a lot of noise in the
% signal)
% % [bg_gxp,std_gxp] = getBackgroundMean(gxp);
% % [bg_gxp2,std_gxp2] = getBackgroundMedian(gxp);
% % thresh_gxp = min(bg_gxp+edgefac*std_gxp,bg_gxp2+edgefac*std_gxp2);
% thresh_gxp = specs.edge_threshold*65535/double(ref_int);
% gxp(gxp<thresh_gxp) = 0;
% % [bg_gxn,std_gxn] = getBackgroundMean(gxn);
% % [bg_gxn2,std_gxn2] = getBackgroundMedian(gxn);
% % thresh_gxn = min(bg_gxn+edgefac*std_gxn,bg_gxn2+edgefac*std_gxn2);
% thresh_gxn = specs.edge_threshold*65535/double(ref_int);
% gxn(gxn<thresh_gxn) = 0;
% % [bg_gyp,std_gyp] = getBackgroundMean(gyp);
% % [bg_gyp2,std_gyp2] = getBackgroundMedian(gyp);
% % thresh_gyp = min(bg_gyp+edgefac*std_gyp,bg_gyp2+edgefac*std_gyp2);
% thresh_gyp = specs.edge_threshold*65535/double(ref_int);
% gyp(gyp<thresh_gyp) = 0;
% % [bg_gyn,std_gyn] = getBackgroundMean(gyn);
% % [bg_gyn2,std_gyn2] = getBackgroundMedian(gyn);
% % thresh_gyn = min(bg_gyn+edgefac*std_gyn,bg_gyn2+edgefac*std_gyn2);
% thresh_gyn = specs.edge_threshold*65535/double(ref_int);
% gyn(gyn<thresh_gyn) = 0;

thresh_edge_min = specs.edge_threshold*65535/double(ref_int);
sumsqr_g = (gxp+gxn).^2+(gyp+gyn).^2;
[grad_cdf,grad_bin] = histcounts(nonzeros(sumsqr_g),nbin,'Normalization','cdf');
thresh_edge_alt = grad_bin(find(grad_cdf>=grad_cdf(1)+0.95*(grad_cdf(end)-grad_cdf(1)),1,'first'));

thresh_edge = max(thresh_edge_min,thresh_edge_alt);

bw_g = imbinarize(sumsqr_g,thresh_edge.^2);
gxp(bw_g==0) = 0;
gxn(bw_g==0) = 0;
gyp(bw_g==0) = 0;
gyn(bw_g==0) = 0;

[bw_edgex,bw_edgey] = fillByEdges(-gxn,gxp,-gyn,gyp,bw_min);
% get AND of two binary images
% get rid of long, thin lines that appear in bw_edge images sometimes
% bw_edgex = bwmorph(bw_edgex,'close',1);
% bw_edgey = bwmorph(bw_edgey,'close',1);
bw_edgex = bwmorph(bw_edgex,'open',1);
bw_edgey = bwmorph(bw_edgey,'open',1);

bw_edge = ones(size(bw_edgex));
bw_edge(bw_edgex==0)=0;
bw_edge(bw_edgey==0)=0;
% fill in gaps in bw_edge by dilating and thining
% the smallest regions will be thrown out in bw_candidates, so it would be
% good to merge small regions that are close together
bw_edge = bwmorph(bw_edge,'clean');
% bw_edge = bwmorph(bw_edge,'dilate',1);
% bw_edge = bwmorph(bw_edge,'thin',1);

% I want to isolate edges that occured in areas with tumors
bw_tumor_edges = bw_edge;
bw_tumor_edges(past_tumors==0) = 0;
bw_tumor_edges(bw_min==0) = 0;
bw_tumor_edges = bwmorph(bw_tumor_edges,'thin',1);

% can in the case of bright tumors, bw_g represents an outline of the
% tumors
bw_edge2 = bw_g;
bw_edge2 = bwmorph(bw_edge2,'clean');
bw_edge2 = imfill(bw_edge2,'holes');
if(ref_int>4e4)
    bw_edge(bw_edge2==1) = 1;
end

% get candidate clusters from bw image
bw_candidates = bw_brightest;
bw_candidates(bw_edge==1)=1;
bw_candidates(bw_signal==1)=1;
bw_candidates(bw_otsu==1)=1;

% clean up candidate 
% bw_candidates = bwmorph(bw_candidates,'dilate',1);
% bw_candidates = bwmorph(bw_candidates,'thin',1);
% bw_candidates = bwmorph(bw_candidates,'erode',1);
% bw_candidates = bwmorph(bw_candidates,'dilate',1);
bw_candidates = bwareaopen(bw_candidates,3);

% make candidate_regions mask by dilating bw_candidates
candidate_regions = bwlabel(bwmorph(bw_candidates,'dilate',specs.clustering_threshold));
candidate_regions = imfill(candidate_regions,'holes');
nreg = max(candidate_regions(:));
candidates = regionprops(logical(candidate_regions));
% calculate stats for candidates and create new bw image
% bw_thresh is a locally thresholded image so that hopefully, we can still
% pick out metastases that are less bright than the primary
bw_thresh = false(size(fp));
thresh_lowest = specs.bw_threshold_low*65535/double(ref_int);
for nr = 1:nreg
    region = fp;
    region(candidate_regions~=nr)=0;
    region_list = nonzeros(region);
    bright_list = nonzeros(region(bw_brightest>0));
    if(isempty(bright_list))
        bright_list = max(region(:));
    end
    dim_list = nonzeros(region(bw_candidates-bw_brightest>0));
    dark_list = nonzeros(region((bw_candidates+bw_brightest+bw_edge)==0));
    % 1) check that there is a bright region to detect
    bri_avg = mean(bright_list);
    %     bri_std = std(bright_list);
    if(~isempty(dark_list))
        drk_avg = mean(dark_list);
        drk_std = std(dark_list);
    else
        drk_avg = 0;
        drk_std = 0;
    end
    % need to decide is this region has a true signal for tumor
    % a) the dark region needs to be distinct from the bright
    candidates(nr).isSeparable = bri_avg > drk_avg+2*drk_std;
    % b) the region should have an edge, which indicates a sharply
    % delineated bright spot in the region, or one of the brightest pixels,
    % - otherwise, it might realy just be all noise
    candidates(nr).edgePixelCount = sum(bw_edge(region~=0));
    candidates(nr).brightPixelCount = sum(bw_brightest(region~=0));
    candidates(nr).isSaturated = false;
    % 2) obtain threshold to separate bright and dark region if this is a
    % good candidate - otherwise
    if(candidates(nr).isSeparable && ...
        candidates(nr).edgePixelCount+candidates(nr).brightPixelCount > 0)
        % we want to keep most of the pixels in the bright list - however,
        % some will be unusually low because bw_edge can include dark
        % pixels adjacent to bright ones - that is why I am using the cdf
        % to pick a threshold for the upper 95% of pixels; however, when
        % there is only a small number of pixels, this method ends up
        % including all the pixels, even the dimmest. To counter this, I am
        % removing pixels dimmer than the average dark pixel + 2*std
        
        % check for saturation
        if(length(bright_list)>1000 && double(ref_int)>1000)
            [bri_cdf,bri_bin] = histcounts(bright_list,1000,'Normalization','cdf');
            per_gt95 = 1-bri_cdf(find(bri_bin>0.95,1));
            if(per_gt95>0.4)
                candidates(nr).isSaturated = true;
            end
        end
        
        if(candidates(nr).isSaturated)
            % Otsu's method should be pretty reasonable at separating
            % saturated signal from unsaturated everything else
            candidates(nr).threshold = 0.9*graythresh(region_list);
        else
            % another threshold method is to use Otsu's method on the dim_list
            % since this should contain some dark and some bright pixels
            if(length(dim_list)>1)
                % Otsu's method needs at least two points to make sense
                thresh_dim = graythresh(dim_list);
                % check to see that this threshold is not too low compared to
                % dark_list
                if(~isempty(dark_list))
                    [drk_cdf,drk_bin] = histcounts(dark_list,1000,'Normalization','cdf');
                    thresh_drk99 = drk_bin(find(drk_cdf>0.99,1));
                    thresh_dim = max(thresh_dim,thresh_drk99);
                end
            else
                thresh_dim = thresh_highest;
            end
            
            thresh_high = max(thresh_signal,thresh_otsu);
            candidates(nr).threshold = min(thresh_dim,thresh_high);
            % low limit to screen out noise
            candidates(nr).threshold = max(thresh_lowest,candidates(nr).threshold);
            % high limit to include bright signals
            candidates(nr).threshold = min(thresh_highest,candidates(nr).threshold);
        end
    else
        if(~isempty(region_list))
            candidates(nr).threshold = max(region_list(:));
        else
            candidates(nr).threshold = 0;
        end
    end
    
    % 3) test sub regions detected using edge detection and add those that
    % pass the test to bw_thresh
    bw_edge_reg = bwmorph(bw_edge,'dilate',2);
    bw_edge_reg(bw_min==0) = 0;
    bw_edge_reg(region==0) = 0;
    bw_edge_reg = bwlabel(bw_edge_reg);
    n_edge = max(bw_edge_reg(:));

    for ie = 1:n_edge
        sub_bw = false(size(bw_edge_reg));
        sub_bw(bw_edge_reg==ie) = 1;
        sub_region = region;
        sub_region(sub_bw==0) = 0;
        edge_reg = sub_region(sub_bw==1);
        thresh_eg = graythresh(edge_reg);
        % will accept as real peak if the bright points form a peak (are
        % centered inside the region) 
        bright_eg = imbinarize(sub_region,thresh_eg);
        dim_eg = sub_bw-bright_eg;
        filled_eg = imfill(dim_eg,'holes');
        if(bwarea(filled_eg)>=bwarea(sub_bw))
            % hole was filled and peak test was passed
            % for large tumors, want to keep more stuff
            if(length(edge_reg)>100)
                bw_thresh(bwmorph(sub_bw,'thin',2)==1) = 1;
            else
                bw_thresh(bright_eg==1) = 1;
            end
        end
        %         imshow(sub_region);
        %         plotBWOutline_v1_0(bright_eg);
        %         plotBWOutline_v1_0(sub_bw);
    end
    
    % 4) create BW image made using local thresholds
    bw_thresh(region>candidates(nr).threshold) = 1;
    % add edges back in since some dim tumors might not meet thresholds
    bw_thresh(bw_tumor_edges) = 1;
    candidates(nr).thesholdedArea = sum(region(:)>candidates(nr).threshold);
end
% clean up bw_thresh by first merging close regions using dilate and thin
% and erode to remove single pixel thin regions, then dilate back to size
bw_thresh = bwareaopen(bw_thresh,3);
% bw_thresh = bwmorph(bw_thresh,'dilate',1);
% bw_thresh = bwmorph(bw_thresh,'thin',1);
% bw_thresh = bwmorph(bw_thresh,'erode',1);
% bw_thresh = bwmorph(bw_thresh,'dilate',1);

% any pixel on bright list and between two edges will be assumed to be
% signal rather than noise
bw_conservative = ones(size(fp));
bw_conservative(bw_brightest==0)=0;
bw_conservative(bw_edge==0)=0;

% merge bw_thresh and bw_conservative
bw = bw_conservative;
bw(bw_thresh==1)=1;
bw = bwareaopen(bw,3);

% dilate to help merge close structures and thin to return to size;
bw = bwmorph(bw,'dilate',0.5*specs.clustering_threshold);
bw = bwmorph(bw,'thin',0.5*specs.clustering_threshold);
% fill holes
bw = imfill(bw,'holes');

% remove 1's inside eye and outside body
bw(eye==1)=0;
bw(body==0)=0;
end

function [avg_g,std_g] = getBackgroundMean(g)
g_list = nonzeros(g);
avg1 = mean(g_list);
std1 = std(g_list);

g_list2 = g_list(g_list<avg1+2*std1);
avg_g = mean(g_list2);
std_g = std(g_list2);
end

function [med_g,std_g] = getBackgroundMedian(g)
g_list = nonzeros(g);
[cdf,bin] = histcounts(g_list,1000,'Normalization','cdf');
med_g = bin(find(cdf>=0.5,1));
std_est1 = bin(find(cdf>=0.8413,1))-med_g;
std_est2 = (bin(find(cdf>=0.9772,1))-med_g)*0.5;
std_g = (std_est1*2+std_est2)/3;
end

function [med_g,std_g] = getZeroMeanBackgroundStd(g)
g_list = nonzeros(g);
[cdf,bin] = histcounts(g_list,1000,'Normalization','cdf');
med_g = 0;
std_est1 = bin(find(cdf>=0.3413,1))-med_g;
std_est2 = (bin(find(cdf>=0.4772,1))-med_g)*0.5;
std_est3 = (bin(find(cdf>=0.4986,1))-med_g)/3;
std_g = (std_est1+std_est2+std_est3)/3;
end

function [avg_g,std_g] = getBackgroundMeanOrZero(g)
g_list = nonzeros(g);
avg1 = mean(g_list);
std1 = std(g_list);

g_list2 = g_list(g_list<avg1+2*std1);

% if the distribution is as expected, the vast majority of the values are
% part of a normal distribution and there is a long tail on the high end
% Since g_list2 contains values lower than 2 std deviations from the mean,
% it should contain 97.7% of the values (95.4/2+50)
% So a test that we can use to see if the distribution is what we expect is
% to make sure g_list2 contains higher than a certain percentage of the
% input values
if(length(g_list2)/length(g_list) > 0.9)
    avg_g = mean(g_list2);
    std_g = std(g_list2);
else
    % looks like there isn't a discrete destribution of background noise
    % will assume everything is signal
    avg_g = 0;
    std_g = 0;
end
end

function [marked] = markValleyToPeak(valleys,peaks,line_min)
% input is a non-positive curve (valleys) and a non-negative curve (peaks)
% The idea is to find the extrema (valleys and peaks) of section of curve 
% (a section is surrounded be zeros on either end). marked will be a binary
% line where regions with a valley to the left and a peak to the right are
% marked with ones and everything else is a zero

% max peak to valley distance for pairing off peaks and valleys
max_pvdist = 20;

% 1) find valleys and peaks
valley_starts = find((valleys(1:end-1)==0).*((valleys(2:end)-valleys(1:end-1))~=0));
valley_locs = zeros(1,length(valley_starts));
for i = 1:length(valley_starts)-1
    [~,loc] = min(valleys(valley_starts(i):valley_starts(i+1)-1));
    valley_locs(i) = loc-1+valley_starts(i);
end
% do final iteration separately since range is to end of line
vlen = length(valley_starts);
if(vlen > 0)
    [~,loc] = min(valleys(valley_starts(vlen):end));
    valley_locs(vlen) = loc-1+valley_starts(vlen);
end

peak_starts = find((peaks(1:end-1)==0).*((peaks(2:end)-peaks(1:end-1))~=0));
peak_locs = zeros(1,length(peak_starts));
for i = 1:length(peak_starts)-1
    [~,loc] = max(peaks(peak_starts(i):peak_starts(i+1)-1));
    peak_locs(i) = loc-1+peak_starts(i);
end
% do final iteration separately since range is to end of line
plen = length(peak_starts);
if(plen > 0)
    [~,loc] = max(peaks(peak_starts(plen):end));
    peak_locs(plen) = loc-1+peak_starts(plen);
end

% 2) pair up valleys and peaks to make marked
marked = zeros(1,length(valleys));
if(vlen > 0)
    next_valley = 1;
    next_peak = find(peak_locs>valley_locs(next_valley),1);
    while(next_valley < vlen && ~isempty(next_peak))
        if(valley_locs(next_valley+1) >= peak_locs(next_peak) && ...
                peak_locs(next_peak)-valley_locs(next_valley) <= max_pvdist)
            % case 1: the closest extrema after the current valley is a
            % peak that is within max_pvdist
            marked(valley_locs(next_valley):peak_locs(next_peak)) = ...
                line_min(valley_locs(next_valley):peak_locs(next_peak));
        end
        % case 2: the closest extrema after the current valley is a valley and no marking needed
        % iterate
        next_valley = next_valley+1;
        next_peak = find(peak_locs>valley_locs(next_valley),1);
    end
    % handle final valley as long as vlen > 0 and there is another peak
    if(~isempty(next_peak))
        % variables have already been iterated
        marked(valley_locs(next_valley):peak_locs(next_peak)) = ...
            line_min(valley_locs(next_valley):peak_locs(next_peak));
    end
end
    
end

function [bw] = fillHorizontalEdges(xneg,xpos,bw_min)
% xneg is non-positive and xpos is non-negative
    bw = zeros(size(xneg));
    for row = 1:size(bw,1)
        bw(row,:) = markValleyToPeak(xneg(row,:),xpos(row,:),bw_min(row,:));
    end
end

function [bw_x,bw_y] = fillByEdges(xneg,xpos,yneg,ypos,bw_min)
% xneg is non-positive and xpos is non-negative and similar for yneg and
% ypos
    bw_x = fillHorizontalEdges(xneg,xpos,bw_min);
    bw_y = fillHorizontalEdges(yneg',ypos',bw_min')';
end
