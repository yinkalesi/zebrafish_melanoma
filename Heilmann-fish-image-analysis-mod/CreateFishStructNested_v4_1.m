%% CreateFishStructFunction
%  Version 4.0
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 4/2/18
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Version History
%  4.0: updated some of the models based on revisions on  Silja data
%  *    modified CreateFishStruct_v4_0 to a function to run in parfor loop
%  4.1: allows only analyzing a subset of the fish

% ***MODIFIED using CreateNestedCFS_v1_0***
function [tumors,totals] = CreateFishStructNested_v4_1(start_step,fish_to_skip,show,batch,baseDir,timePoint,fish_specs,...
img_names,gender_file,lookup_file,elimination_file,implantSize,implantSizeName,location,fish_struct_save_name)

%% Start log
diary([fish_struct_save_name '_log.txt']);

%% Set up names, directories and descriptions
start_time = clock; % start timer

% fish_dir1 = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\White Lab\td013017 Adult Transplant ZMEL1 1DPT/NO RAD';
% fish_dir2 = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\White Lab\td013017 Adult Transplant ZMEL1 3DPT/NO RAD';
% fish_dir3 = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\White Lab\td013017 Adult Transplant ZMEL1 5DPT/NO RAD';
% fish_dir4 = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\White Lab\td013017 Adult Transplant ZMEL1 7DPT/NO RAD';
% fish_dir5 = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\White Lab\td013017 Adult Transplant ZMEL1 9DPT/NO RAD';

% fish_struct_save_name = 'FISH_NORAD_031017_struct.mat';
% baseDir = {fish_dir1,fish_dir2,fish_dir3,fish_dir4,fish_dir5}; %names of time point folders
% timePoint = {'1','3','5','7','9'}; 
% GFP_exposure_times = [4000 3000 3000 206.4 300];
% batch    = {'B013017'};% 'BATCH2'}; %names of batch folders
% implantSize = {'5x10^5'};% '5x10^5' '1x10^5'}; % first half of group name folder names
% implantSizeName = {'1e5'};% '5e5' '1e6'};
% location = {'VENTRAL'};% 'DORSAL'}; % second half of group name folder names

% fish_specs = struct('minimum_area',50000,'maximum_length',1350,...
%     'maximum_width',600,'nose_coordinate',[10 300],'eye_diameter',80,...
%     'eye_coordinate',[100 290],'clustering_threshold',10,...
%     'gfp_exposures',B0213_GFP_times,'imageConversionFactor',257,...
%     'eyeBrightnessFactor',20,'rfp_gfp_ratio',1,'minimum_eye_area',1500,...
%     'mean_intensities',[6.22e4,2.00e3,2.00e3],'bw_threshold_high',0.01,...
%     'bw_threshold_low',0.005);

% % setting to display results of image analysis
% show = 0;
% start_step = 5;
% fish_to_skip = [];

% start_step = 1;
step_num = 0;
%% Step 1: load original images in to data structure

step_num = step_num+1;
if(start_step <= step_num)
BYTIME = CreateFishStruct_step1_v2_1(baseDir,timePoint,batch,...
implantSizeName,location,img_names,fish_specs);
end

%% Step 2: rescale image intensities

step_num = step_num+1;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step2_v1_6(fish_specs,show);
end

%% Step 3: flip images

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step3_v1_1(show);
end

%% Step 4: Rough body outline

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step4_v1_3(fish_specs,show);
end

%% Step 5: Rotation and cropping

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step5_v1_2(fish_specs,show);
end

%% Step 6: Better body outline

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step6_v1_4(fish_specs,show);
end

%% Step 7: Eye outline (used to be step 8 so the function called below refers to step 8)

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step8_v1_3(fish_specs,show);
end

%% Step 8: Spine outline (now uses eye position - thus 'step 8' will be done first)

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step7_v2_2(fish_struct_save_name,show);
end

%% Step 9: Left-right transforms

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step9_v2_1(fish_specs,show);
end

%% Step 10: Merge left and right images

step_num = step_num+1;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step10_v1_3();
end

%% Step 10.5: Compare fish and identify over time points

step_num = step_num+0.5;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_identify_v2_1(lookup_file,show);
end

%% Step 10.6: Assign fish genders
step_num = step_num+0.1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_gender_v1_2(gender_file,fish_specs,show);
end

%% Step 11: Figure out how to transform images at a different times

step_num = step_num+0.4;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step11_v2_2(fish_specs,show);
end

%% Step 12: Save transformed image

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step12_v2_1(fish_specs,show);
end

%% Step 13: Substract autofluorescence

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step13_v2_9(fish_specs,fish_struct_save_name,show);
end

%% Step 14: Remove object no longer there at next time point
%  (may be detrimental to my data)

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
% CreateFishStruct_step14_v1_1(show);
end

%% Step 15: Calculate pairwise distances and use to cluster objects

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step15_v1_1(fish_specs,show);
end

%% Step 15.5: Remove fish we don't want to analyze
skipList = fish_to_skip;
for sli = 1:length(skipList)
dati = find(BYTIME(1).dataToUse==skipList(sli),1);
if(~isempty(dati))
BYTIME(1).dataToUse = BYTIME(1).dataToUse([1:dati-1 dati+1:end]);
end

end

%% Step 16: Separating merged objects

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_step16_v1_5(fish_specs,fish_struct_save_name,show);
end

%% Step 16.5: Elimating certain previously selected fish - screened out with there images

step_num = step_num+0.5;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 CreateFishStruct_eliminations_v1_1(elimination_file,show);
end

%% Step 17: Transform to stardard fish shape
%  (not needed for me)

step_num = step_num+0.5;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
%      [ aveFish] = CreateFishStruct_step17_v1_1(fish_specs,show);
end

%% Step 18: Applying transform functions found in previous step to all images
%  (not needed for me)

step_num = step_num+1;
close;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
% CreateFishStruct_step18_v1_1(fish_specs,aveFish,show);
end

%% Step 19: get and save tumor summary data (formerly step 20)

% save the most salient info (the produced structures are also saved in the function)
step_num = step_num+1;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
 [tumors,totals,~] = access_fish_struct_v5_1([fish_struct_save_name '_summary']);
end

%% Step 20: Save data struct (formerly step 19)

step_num = step_num+1;
BYTIME(1).specs = fish_specs;
if(start_step <= step_num)
% ***MODIFIED using CreateNestedCFS_v1_0***
% CreateFishStruct_step19_v1_2(fish_struct_save_name);
end



%% Wrap up
duration = etime(clock,start_time);

durahrs = floor(duration/3600);
duramins = floor((duration - durahrs*3600)/60);
durasecs = round(duration - durahrs*3600 - duramins*60);

fprintf('CreateFishStruct took %i hours, %i minutes and %i seconds\n',...
durahrs,duramins,durasecs);

%% End log
diary off;



% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step2_v1_6(specs,show)
%% CreateFishStruct_step2
%  Version 1.6
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 2/10/20
%  Project: Tumor Growth, Logarithmic Continuum Form

%%  Rescale overall intensity of images images based on cut out of background (BG)
%   and save intensity rescaled images in rounds(1)
%   Note: if all images are taken with the exact same settings and ligthing
%   this step should not be nessesary and may be skipped!
%   1.1: changing rescaling method to increase contrast
%   1.2: not satisfied with 1.1 - 
%   1.3: a) use real exposure times to normalize image brightness, b)
%   substract out average background signal
%   1.4: using bf(:,:,1)/able to handle rgb bf object
%   1.6: minor adjustments to code

tic;

disp('Step 2: Rescaling brightness of BF, GFP and RFP images based on average background intensity...');

rect1 = [400 1 800 110]; % [xmin ymax width height] of background cutout
rect2 = [400 1 800 110];
% rect2 = [188 501 800 110];

% Target background intensities
exposure_ratios = mean(specs.experimental_exposure_times,2)./specs.experimental_exposure_times;

exp_corr = specs.normal_exposure_times./specs.experimental_exposure_times;

% conversion factor from uint8 to uint16, if necessary
uintconv = specs.imageConversionFactor;
meanIntensities = zeros(4,length(BYTIME));

for tt = 1:length(BYTIME)
FISH = BYTIME(tt).fishes;
originalInts = zeros(length(BYTIME(tt).fishes),8);
meanInts = zeros(1,5);
for ff = BYTIME(end).dataAt{tt}
LR = 0;
while LR<2
LR = LR+1;
switch LR
% retrieave images from data struct.
case 1                    
bf  = uint16(FISH(ff).bf1)*uintconv;
gfp = uint16(FISH(ff).gfp1)*uintconv;
rfp = uint16(FISH(ff).rfp1)*uintconv;
rgb = uint16(FISH(ff).rgb1)*uintconv;
TITLE = 'Right side';
rect = rect1;
case 2
bf  = uint16(FISH(ff).bf2)*uintconv;
gfp = uint16(FISH(ff).gfp2)*uintconv;
rfp = uint16(FISH(ff).rfp2)*uintconv;
rgb = uint16(FISH(ff).rgb2)*uintconv;
TITLE = 'Left side';
rect = rect2;
end

% use thresholding to figure out where fish is
bf_thresh = graythresh(bf(:,:,1));
rgb_thresh = graythresh(rgb(:,:,1));
bf_mask = 1-imbinarize(bf(:,:,1),bf_thresh*1.2);
rgb_mask = 1-imbinarize(rgb(:,:,1),rgb_thresh*1.2);
bf_mask = bwmorph(bf_mask,'dilate',20);
rgb_mask = bwmorph(rgb_mask,'dilate',20);
% cut out piece of BG at top edge of image (contains no fish)
if(ismatrix(bf))
cropBF = bf(bf_mask==0);
graybf = bf;
else
graybf = rgb2gray(bf);
cropBF = graybf(bf_mask==0);
end
cropGFP = gfp(bf_mask==0);
cropRFP = rfp(bf_mask==0);
gray = rgb2gray(rgb);
cropRGB = gray(rgb_mask==0);
%             cropBF = imcrop(bf,rect);
%             cropGFP = imcrop(gfp,rect);
%             cropRFP = imcrop(rfp,rect);
%             cropRGB = imcrop(rgb2gray(rgb),rect);

% determine mean and std of BG in BF, GFP and RFP images
MBF    = mean(cropBF(:));
STDBF  = std(double(cropBF(:)));
MGFP   = mean(cropGFP(:));
STDGFP = std(double(cropGFP(:)));
MRFP   = mean(cropRFP(:));
STDRFP = std(double(cropRFP(:)));
MRGB   = mean(cropRGB(:));
STDRGB = std(double(cropRGB(:)));

% get rid of extreme values in the cropped images - example dark
% spot due to random dirt.

temp_cropBF = cropBF(cropBF<MBF+2*STDBF);
temp_cropBF = temp_cropBF(temp_cropBF>MBF-2*STDBF);
temp_cropGFP = cropGFP(cropGFP<MGFP+2*STDGFP);
temp_cropGFP = temp_cropGFP(temp_cropGFP>MGFP-2*STDGFP);
temp_cropRFP = cropRFP(cropRFP<MRFP+2*STDRFP);
temp_cropRFP = temp_cropRFP(temp_cropRFP>MRFP-2*STDRFP);
temp_cropRGB = cropRGB(cropRGB<MRGB+2*STDRGB);
temp_cropRGB = temp_cropRGB(temp_cropRGB>MRGB-2*STDRGB);

% identify which regions were removed in mask
gfp_mask = bf_mask;
rfp_mask = bf_mask;
bf_mask(graybf>=MBF+2*STDBF) = 1;
bf_mask(graybf<=MBF-2*STDBF) = 1;
gfp_mask(gfp>=MGFP+2*STDGFP) = 1;
gfp_mask(gfp<=MGFP-2*STDGFP) = 1;
rfp_mask(rfp>=MRFP+2*STDRFP) = 1;
rfp_mask(rfp<=MRFP-2*STDRFP) = 1;
rgb_mask(gray>=MRGB+2*STDRGB) = 1;
rgb_mask(gray<=MRGB-2*STDRGB) = 1;

% calculate mean BG intensity without extreme values caused by random dirt
MBF = mean(temp_cropBF(:));
MGFP = mean(temp_cropGFP(:));
MRFP = mean(temp_cropRFP(:));
MRGB = mean(temp_cropRGB(:));

% create new image with better contrast
bf_new = uint16(double(bf)*exp_corr(1,tt));
gfp_new = uint16((double(gfp)-MGFP)*exp_corr(2,tt));
rfp_new = uint16((double(rfp)-MRFP)*exp_corr(3,tt));
rgb_new = uint16(double(rgb)*exp_corr(4,tt));
%             disp(['GFP: ' num2str(MGFP) ', RFP: ' num2str(MRFP) ', RFP/GFP: ' num2str(MRFP/MGFP)]);
meanInts = meanInts + [MBF MGFP MRFP MRGB 1];
originalInts(ff,(1:4)+(LR-1)*4) = [MBF MGFP MRFP MRGB];

% save rescaled images in rounds(1)
switch LR
case 1
FISH(ff).rounds(1).bf1  = uint16(bf_new);
FISH(ff).rounds(1).gfp1 = uint16(gfp_new);
FISH(ff).rounds(1).rfp1 = uint16(rfp_new);
FISH(ff).rounds(1).rgb1 = uint16(rgb_new);
case 2
FISH(ff).rounds(1).bf2  = uint16(bf_new);
FISH(ff).rounds(1).gfp2 = uint16(gfp_new);
FISH(ff).rounds(1).rfp2 = uint16(rfp_new);
FISH(ff).rounds(1).rgb2 = uint16(rgb_new);
end

if show==1
%visually inspect intensity normalisation/rescaling of brigthness
figure(1);
%                 h = figure(1);
%                 set(h,'Position',[10,50,1900,950]);
subplot(2,2,1)
imshowpair(bf_new,1-bf_mask);
title(['Rescaled BF. Fish number ' num2str(ff) ',  time: ' num2str(tt)],'FontSize',13)
%                 % plot rect box;
%                 hold on;
%                 plot([rect1(1) rect1(1)+rect1(3) rect1(1)+rect1(3) rect1(1) rect1(1)],...
%                      [rect1(2) rect1(2) rect1(2)+rect1(4) rect1(2)+rect1(4) rect1(2)]);
%                 plot([rect2(1) rect2(1)+rect2(3) rect2(1)+rect2(3) rect2(1) rect2(1)],...
%                      [rect2(2) rect2(2) rect2(2)+rect2(4) rect2(2)+rect2(4) rect2(2)]);
%                 hold off;
outline1 = bwareaopen(gfp_mask,9);
subplot(2,2,2)
maxval = max([gfp_new(:); rfp_new(:)]);
imshow(gfp_new,[0 maxval]);
plotBWOutline_v1_0(outline1,'y-');
title(['Rescaled GFP. Fish number ' num2str(ff) ',  time: ' num2str(tt)],'FontSize',13)

outline2 = bwareaopen(rfp_mask,9);
subplot(2,2,3)
imshow(rfp_new,[0 maxval]);
plotBWOutline_v1_0(outline2,'y-');
title(['Rescaled RFP. Fish number ' num2str(ff) ',  time: ' num2str(tt)],'FontSize',13)

subplot(2,2,4)
imshowpair(rgb_new,1-rgb_mask);
title(['Rescaled RGB. Fish number ' num2str(ff) ',  time: ' num2str(tt)],'FontSize',13)

pause(0.1);
end
end
end
BYTIME(tt).fishes = FISH;
BYTIME(tt).originalIntensities = originalInts;
meanInts = meanInts/meanInts(end);
BYTIME(tt).meanIntensities = meanInts(1:4);
meanIntensities(:,tt) = meanInts(1:4)';
disp(['Average Intensities (BF, GFP, RFP, Color): ' num2str(meanInts(1:4))]);
end
BYTIME(1).meanIntensitySummary = meanIntensities;
disp('Exposure ratios: ');
fprintf(' BF_er  BF_mi GFP_er GFP_mi RFP_er RFP_mi Col_mi Col_er \n');
trange = 1:length(BYTIME);
fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n',...
[exposure_ratios(1,trange); mean(meanIntensities(1,:))./meanIntensities(1,:);...
exposure_ratios(2,trange); mean(meanIntensities(2,:))./meanIntensities(2,:);...
exposure_ratios(3,trange); mean(meanIntensities(3,:))./meanIntensities(3,:);...
exposure_ratios(4,trange); mean(meanIntensities(4,:))./meanIntensities(4,:)]);

toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step3_v1_1(show)
%% CreateFishStruct_step3
%  Version 1.0
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 3/10/16
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Flip images of fish which has nose to the right and save in rounds(1)

%% Version History
%  1.1: works with rgb images

disp('Step 3: Flipping images of fish which has nose to the right and save in rounds(1)...');

tic;

for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}

% retrieave images from data struct.
gfp1 = FISH(ff).rounds(1).gfp1;
rfp1 = FISH(ff).rounds(1).rfp1;
bf1  = FISH(ff).rounds(1).bf1;
rgb1 = FISH(ff).rounds(1).rgb1;
bf2  = FISH(ff).rounds(1).bf2;

% flip red green and blue channels of rgb1 image one at a time
temp1 = rgb1;
temp1(:,:,1) = fliplr(rgb1(:,:,1));
temp1(:,:,2) = fliplr(rgb1(:,:,2));
temp1(:,:,3) = fliplr(rgb1(:,:,3));

% save in rounds(1)
FISH(ff).rounds(1).gfp1 = fliplr(gfp1);
FISH(ff).rounds(1).rfp1 = fliplr(rfp1);
FISH(ff).rounds(1).bf1  = fliplr(bf1);
FISH(ff).rounds(1).rgb1 = temp1;

if show==1
%visually inspect that all fish face in right direction!
figure(1)
imshowpair(imresize(bf2,0.5),imresize(FISH(ff).rounds(1).bf1,0.5))
hold on
title(['Flipping images where fish nose where facing right. Green right side, pink left side. (Both noses should be facing left!). Fish number: ' num2str(ff) ',  time: ' num2str(t)],'FontSize',11)
pause(1);
end
end
BYTIME(t).fishes = FISH;
end

toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step4_v1_3(specs,show)
%% CreateFishStruct_step4
%  Version 1.3
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 3/13/18
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Find approx. outline of fish bodies
% this rough outline/boundary of fish body will be used later for rotating the
% fish, translating the 'nose' to a specific commen point and for cropping
% the images, and for finding a more accurate outline later. 
%% Version History:
%  1.1: changing how tail and fin bones are removed from body outline
%  1.2: processing rgb images
%  1.3: adjustments for Silja data

disp('Step 4: Finding approx outline of fish bodies and save in rounds(1)...');


tic ;

% 50000 is a reasonable number
min_conn_region_size = specs.minimum_area; % this number needs to be slightly less than the size of the smallest fish in the data set


ww = specs.maximum_length; % slightly longer than length of longest fish
hh = specs.maximum_width; % slightly wider than width of fattest/widest (or most not horizontal) fish

for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}

% retrieave images from data struct.
bf1  = FISH(ff).rounds(1).bf1;
bf2  = FISH(ff).rounds(1).bf2;
rgb1  = FISH(ff).rounds(1).rgb1;
rgb2  = FISH(ff).rounds(1).rgb2;

% merge rgb channels into on grayscale image
%         
% for color images, use green field:
if(ndims(bf1)==3)
bw1 = 1-imbinarize(bf1(:,:,2));
bw2 = 1-imbinarize(bf2(:,:,2));
else
%             thresh1 = max(mean(double(bf1(:))/65535),graythresh(bf1));
thresh1 = graythresh(bf1);
bw1 = 1-imbinarize(bf1,thresh1);
%             thresh2 = max(mean(double(bf2(:))/65535),graythresh(bf2));
thresh2 = graythresh(bf2);
bw2 = 1-imbinarize(bf2,thresh2);
end
cbw1 = 1-imbinarize(rgb1(:,:,2));
cbw2 = 1-imbinarize(rgb2(:,:,2));

% some of the bones in the fins show up as thin lines - will try to
% remove them by eroding and dilating to cutoff thin regions
bw1 = imfill(bw1,'holes');
bw1 = bwmorph(bw1,'thin',10);
bw1 = bwmorph(bw1,'erode',1);
bw1 = bwmorph(bw1,'dilate',1);
bw1 = bwmorph(bw1,'thicken',10);
bw2 = imfill(bw2,'holes');
bw2 = bwmorph(bw2,'thin',10);
bw2 = bwmorph(bw2,'erode',1);
bw2 = bwmorph(bw2,'dilate',1);
bw2 = bwmorph(bw2,'thicken',10);
cbw1 = bwmorph(cbw1,'thin',10);
cbw1 = bwmorph(cbw1,'erode',1);
cbw1 = bwmorph(cbw1,'dilate',1);
cbw1 = bwmorph(cbw1,'thicken',10);
cbw2 = bwmorph(cbw2,'thin',10);
cbw2 = bwmorph(cbw2,'erode',1);
cbw2 = bwmorph(cbw2,'dilate',1);
cbw2 = bwmorph(cbw2,'thicken',10);
% tail needs more work (see step 6)

% fill in holes in connected regions of 1's
bw1 =imfill(bw1,'holes');
bw2 =imfill(bw2,'holes');
cbw1 =imfill(cbw1,'holes');
cbw2 =imfill(cbw2,'holes');

% remove connected regions of 1's with fewer than min_size pixels
bw1 = bwareaopen(bw1,min_conn_region_size);
bw2 = bwareaopen(bw2,min_conn_region_size);
cbw1 = bwareaopen(cbw1,min_conn_region_size);
cbw2 = bwareaopen(cbw2,min_conn_region_size);

% find coordinates of fish 'nose' on left and right side images
[y1,x1] = find(bw1,5,'first');
noseCoor1 = [round(mean(x1)) round(mean(y1))];
[y2,x2] = find(bw2,5,'first');
noseCoor2 = [round(mean(x2)) round(mean(y2))];
[y3,x3] = find(cbw1,5,'first');
noseCoor3 = [round(mean(x3)) round(mean(y3))];
[y4,x4] = find(cbw2,5,'first');
noseCoor4 = [round(mean(x4)) round(mean(y4))];

% remove 1's in regions beyond the size of mask
bw1(1:noseCoor1(2)-hh,:) = 0;
bw2(1:noseCoor2(2)-hh,:) = 0;
bw1(noseCoor1(2)+hh:end,:) = 0;
bw2(noseCoor2(2)+hh:end,:) = 0;
bw1(:,1:noseCoor1(1)-10,:) = 0;
bw2(:,1:noseCoor2(1)-10,:) = 0;
bw1(:,noseCoor1(1)+ww:end) = 0;
bw2(:,noseCoor2(1)+ww:end) = 0;
cbw1(1:noseCoor3(2)-hh,:) = 0;
cbw2(1:noseCoor4(2)-hh,:) = 0;
cbw1(noseCoor3(2)+hh:end,:) = 0;
cbw2(noseCoor4(2)+hh:end,:) = 0;
cbw1(:,1:noseCoor3(1)-10,:) = 0;
cbw2(:,1:noseCoor4(1)-10,:) = 0;
cbw1(:,noseCoor3(1)+ww:end) = 0;
cbw2(:,noseCoor4(1)+ww:end) = 0;

% make boundary of fish body smoother
bw1 = smoothBW_mod(bw1,10);
bw2 = smoothBW_mod(bw2,10);  
cbw1 = smoothBW_mod(cbw1,10);
cbw2 = smoothBW_mod(cbw2,10); 

if show==1
% visually inspect that approx. fish outlines look ok
figure(1)
subplot(2,1,1)
imshowpair(imresize(bw1,0.5),imresize(bw2,0.5));
hold on
plot(noseCoor1(1).*0.5,noseCoor1(2).*0.5,'*','MarkerSize',15);
plot(noseCoor2(1).*0.5,noseCoor2(2).*0.5,'*','MarkerSize',15);
title(['Finding approx outline of fish bodies. Green right side, pink left side.  Fish number: ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
subplot(2,1,2)
imshowpair(imresize(cbw1,0.5),imresize(cbw2,0.5));
hold on
plot(noseCoor3(1).*0.5,noseCoor3(2).*0.5,'*','MarkerSize',15);
plot(noseCoor4(1).*0.5,noseCoor4(2).*0.5,'*','MarkerSize',15);
pause(0.1);
end

% save approx body mask in rounds(1)
FISH(ff).rounds(1).body1 = logical(bw1);
FISH(ff).rounds(1).body2 = logical(bw2);
FISH(ff).rounds(1).body3 = logical(cbw1);
FISH(ff).rounds(1).body4 = logical(cbw2);
if(sum(FISH(ff).rounds(1).body1(:))==0 || sum(FISH(ff).rounds(1).body2(:))==0)
error('No body outline; need to adjust min_conn_region_size.');
end
end
BYTIME(t).fishes = FISH;
end
toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step5_v1_2(specs,show)
%% CreateFishStruct_step5
%  Version 1.2
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 10/18/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Rotate fish body to horizontal, trabslate such that nose is at fixed coor and crop images
% using rough fish body outline found above
% save rotated and cropped images in rounds(2)
%% Version History
%  9/5/17: changed width and height of rect to ww-1 and hh-1, since the
%  dimension of cropped image is the number + 1
%  1.1: switched to find_rot_angle_mod2, changed imrotate parameter to
%  'loose' so no part of fish is accidentally cropped, saving rotation
%  angle for future reference;
%  1.2: processing rgb image

disp('Step 5: Rotating translating and cropping all images and save in rounds(2)...');

tic;

close all

ww = specs.maximum_length; % width of cropped images
hh = specs.maximum_width;  % heigth of cropped images
nose_coor_fix = specs.nose_coordinate;

for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}
LR = 0;
while LR<2
LR = LR+1;
switch LR
% retrieave images from data struct.
case 1
bf    = FISH(ff).rounds(1).bf1;
gfp   = FISH(ff).rounds(1).gfp1;
rfp   = FISH(ff).rounds(1).rfp1;
rgb   = FISH(ff).rounds(1).rgb1;
body  = FISH(ff).rounds(1).body1;
cbody = FISH(ff).rounds(1).body3;
case 2
bf    = FISH(ff).rounds(1).bf2;
gfp   = FISH(ff).rounds(1).gfp2;
rfp   = FISH(ff).rounds(1).rfp2;
rgb   = FISH(ff).rounds(1).rgb2;
body  = FISH(ff).rounds(1).body2;
cbody = FISH(ff).rounds(1).body4;
end

% find angle of fish body with respect to horizontal
if(sum(body(:))==0)
error('No fish body detected for f: %i at t: %i',ff,t);
end
ang = find_rot_angle_mod2(body,specs.minimum_area);
ang2 = find_rot_angle_mod2(cbody,specs.minimum_area);

% rotate images so that fish becomes horizontal
bodyN = logical(imrotate(body,ang,'nearest','loose'));
bfN   = imrotate(bf,ang,'bicubic','loose');
% want to make it so that the background upon rotation is white
% instead of black for bf image
% will rotate a rectangle and everyting outide the rectangle in
% the rotated image is background
blank = imrotate(true(size(bf)),ang,'nearest','loose');
bfN(blank==0)=mean(bf(bf>mean(bf(:))));
gfpN  = imrotate(gfp,ang,'bicubic','loose');
rfpN  = imrotate(rfp,ang,'bicubic','loose');
cbodyN = logical(imrotate(cbody,ang2,'nearest','loose'));
rgbN = imrotate(rgb,ang2,'bicubic','loose');
blank2 = imrotate(true(size(rgb(:,:,1))),ang2,'nearest','loose');
temp = rgbN(:,:,1);
temp(blank2==0)=mean(temp(temp>mean(temp(:))));
rgbN(:,:,1) = temp;
temp = rgbN(:,:,2);
temp(blank2==0)=mean(temp(temp>mean(temp(:))));
rgbN(:,:,2) = temp;
temp = rgbN(:,:,3);
temp(blank2==0)=mean(temp(temp>mean(temp(:))));
rgbN(:,:,3) = temp;

% get shifted and cropped images            
% find nose coor
[y,x] = find(bodyN,10,'first');
nose_coor = round([mean(x),mean(y)]);
[y2,x2] = find(cbodyN,10,'first');
nose_coor2 = round([mean(x2),mean(y2)]);

% rectangle for cropping [upper left corner coor , width height] -
% makes sure nose is at (10,nose_coor_fix(2))
han = min(nose_coor_fix(2),hh)-1; %height above nose
rect = [nose_coor(1)-nose_coor_fix(1)+1 nose_coor(2)-han ww-1 hh-1];%
rect2 = [nose_coor2(1)-nose_coor_fix(1)+1 nose_coor2(2)-han ww-1 hh-1];

body_temp = shiftAndCut(bodyN,rect,0);
bf_temp   = shiftAndCut(bfN,rect);
gfp_temp  = shiftAndCut(gfpN,rect);
rfp_temp  = shiftAndCut(rfpN,rect);
cbody_temp = shiftAndCut(cbodyN,rect2);
rgb_temp = shiftAndCut(rgbN,rect2);

%             if(rect(1)<1)
% %                 disp(['Shifting right FISH ' num2str(ff) ' at time ' num2str(t)]);
%                 % need to shift image so it can be cropped properly
%                 shift = ceil(-rect(1))+1;
%                 bodyN(:,shift+1:end) = bodyN(:,1:end-shift);
%                 bodyN(:,1:shift) = 0;
%                 right_edge = bfN(:,end-shift+1:end);
%                 bfN(:,shift+1:end) = bfN(:,1:end-shift);
%                 bfN(:,1:shift) = right_edge;
%                 right_edge = gfpN(:,end-shift+1:end);
%                 gfpN(:,shift+1:end) = gfpN(:,1:end-shift);
%                 gfpN(:,1:shift) = right_edge;
%                 right_edge = rfpN(:,end-shift+1:end);
%                 rfpN(:,shift+1:end) = rfpN(:,1:end-shift);
%                 rfpN(:,1:shift) = right_edge;
%                 rect(1) = rect(1)+shift;
%             end
%             if(rect(1)+rect(3)>size(bodyN,2))
% %                 disp(['Shifting left FISH ' num2str(ff) ' at time ' num2str(t)]);
%                 % need to shift image so that entire crop is within image
%                 shift = ceil(rect(1)+rect(3)-size(bodyN,2));
%                 bodyN(:,1:end-shift) = bodyN(:,1+shift:end);
%                 bodyN(:,end-shift+1:end) = 0;
%                 left_edge = bfN(:,1:shift);
%                 bfN(:,1:end-shift) = bfN(:,1+shift:end);
%                 bfN(:,end-shift+1:end) = left_edge;
%                 left_edge = gfpN(:,1:shift);
%                 gfpN(:,1:end-shift) = gfpN(:,1+shift:end);
%                 gfpN(:,end-shift+1:end) = left_edge;
%                 left_edge = rfpN(:,1:shift);
%                 rfpN(:,1:end-shift) = rfpN(:,1+shift:end);
%                 rfpN(:,end-shift+1:end) = left_edge;
%                 rect(1) = rect(1)-shift;
%             end
%             if(rect(2) < 1)
%                 % shift image to allow appropriate cropping
% %                 disp(['Shifting down FISH ' num2str(ff) ' at time ' num2str(t)]);
%                 shift = ceil(-rect(2))+1;
%                 bottom_edge = bodyN(end-shift+1:end,:);
%                 bodyN(shift+1:end,:) = bodyN(1:end-shift,:);
%                 bodyN(1:shift,:) = bottom_edge;
%                 bottom_edge = bfN(end-shift+1:end,:);
%                 bfN(shift+1:end,:) = bfN(1:end-shift,:);
%                 bfN(1:shift,:) = bottom_edge;
%                 bottom_edge = gfpN(end-shift+1:end,:);
%                 gfpN(shift+1:end,:) = gfpN(1:end-shift,:);
%                 gfpN(1:shift,:) = bottom_edge;
%                 bottom_edge = rfpN(end-shift+1:end,:);
%                 rfpN(shift+1:end,:) = rfpN(1:end-shift,:);
%                 rfpN(1:shift,:) = bottom_edge;
%                 rect(2) = rect(2)+shift;
%             end
%             if(rect(2)+rect(4) > size(bodyN,1))
%                 % shift image
% %                 disp(['Shifting up FISH ' num2str(ff) ' at time ' num2str(t)]);
%                 shift = ceil(rect(2)+rect(4)-size(bodyN,1));
%                 top_edge = bodyN(1:shift,:);
%                 bodyN(1:end-shift,:) = bodyN(shift+1:end,:);
%                 bodyN(end-shift+1:end,:) = top_edge;
%                 top_edge = bfN(1:shift,:);
%                 bfN(1:end-shift,:) = bfN(shift+1:end,:);
%                 bfN(end-shift+1:end,:) = top_edge;
%                 top_edge = gfpN(1:shift,:);
%                 gfpN(1:end-shift,:) = gfpN(shift+1:end,:);
%                 gfpN(end-shift+1:end,:) = top_edge;
%                 top_edge = rfpN(1:shift,:);
%                 rfpN(1:end-shift,:) = rfpN(shift+1:end,:);
%                 rfpN(end-shift+1:end,:) = top_edge;
%                 rect(2) = rect(2)-shift;
%             end
%       
%             % output images (rotated and cropped)
%             body_temp = imcrop(bodyN,rect);
%             bf_temp   = imcrop(bfN,rect);
%             gfp_temp  = imcrop(gfpN,rect);
%             rfp_temp  = imcrop(rfpN,rect);

if(any(size(body_temp)-[hh ww]))
error('Size misassignment at f: %i, t: %i',ff,t);
end

if show == 1
% visually inspect that fish where rotated and images cropped
% correctly:
figure(1)
subplot(2,2,(LR-1)*2+1); imshowpair(body_temp,cbody_temp)
title(['Rotate and crop.  Fish: ' num2str(ff) ',  time: ' num2str(t)],'FontSize',11)
subplot(2,2,(LR-1)*2+2); imshowpair(bf_temp,rgb_temp)
pause(0.1);
end

switch LR
case 1
FISH(ff).rounds(2).bf1   = uint16(bf_temp);
FISH(ff).rounds(2).gfp1  = uint16(gfp_temp);
FISH(ff).rounds(2).rfp1  = uint16(rfp_temp);
FISH(ff).rounds(2).rgb1  = uint16(rgb_temp);
FISH(ff).rounds(2).R_body1 = logical(body_temp); 
FISH(ff).rounds(2).R_body3 = logical(cbody_temp);
FISH(ff).transforms.rotation1 = ang;
FISH(ff).transforms.rotation3 = ang2;
FISH(ff).transforms.rotation1_distortion = sum(body_temp(:))/sum(body(:));
FISH(ff).transforms.rotation3_distortion = sum(cbody_temp(:))/sum(cbody(:));
case 2
FISH(ff).rounds(2).bf2   = uint16(bf_temp);
FISH(ff).rounds(2).gfp2  = uint16(gfp_temp);
FISH(ff).rounds(2).rfp2  = uint16(rfp_temp);
FISH(ff).rounds(2).rgb2  = uint16(rgb_temp);
FISH(ff).rounds(2).R_body2 = logical(body_temp);
FISH(ff).rounds(2).R_body4 = logical(cbody_temp);
FISH(ff).transforms.rotation2 = ang;
FISH(ff).transforms.rotation4 = ang2;
FISH(ff).transforms.rotation2_distortion = sum(body_temp(:))/sum(body(:));
FISH(ff).transforms.rotation4_distortion = sum(cbody_temp(:))/sum(cbody(:));
end
end
end
BYTIME(t).fishes = FISH;
end

toc; 
disp(' ')

end

function [img] = shiftAndCut(img,rect,fill)
% crop image to size; if the cropping rectangle is out of frame, shift the
% image and the cropping rectangle to make it in frame
if(nargin > 2)
no_fill = 0;
else
no_fill = 1;
fill = 0;
end

% check for multiple spectra images
if(ndims(img)==3)
dimnum = size(img,3);
else
dimnum = 1;
end
shifts = zeros(1,4);
for dim = 1:dimnum
im = img(:,:,dim);
if(rect(1)<1)
shift = ceil(-rect(1))+1;
if(no_fill)
right_edge = im(:,end-shift+1:end);
else
right_edge = ones(size(im(:,end-shift+1:end)))*fill;
end
im(:,shift+1:end) = im(:,1:end-shift);
im(:,1:shift) = right_edge;
shifts(1) = shift;
end
if(rect(1)+rect(3)>size(im,2))
shift = ceil(rect(1)+rect(3)-size(im,2));
if(no_fill)
left_edge = im(:,1:shift);
else
left_edge = ones(size(im(:,1:shift)))*fill;
end
im(:,1:end-shift) = im(:,1+shift:end);
im(:,end-shift+1:end) = left_edge;
shifts(2) = shift;
end
if(rect(2) < 1)
shift = ceil(-rect(2))+1;
if(no_fill)
bottom_edge = im(end-shift+1:end,:);
else
bottom_edge = ones(size(im(end-shift+1:end,:)))*fill;
end
im(shift+1:end,:) = im(1:end-shift,:);
im(1:shift,:) = bottom_edge;
shifts(3) = shift;
end
if(rect(2)+rect(4) > size(im,1))
shift = ceil(rect(2)+rect(4)-size(im,1));
if(no_fill)
top_edge = im(1:shift,:);
else
top_edge = ones(size(im(1:shift,:)))*fill;
end
im(1:end-shift,:) = im(shift+1:end,:);
im(end-shift+1:end,:) = top_edge;
shifts(4) = shift;
end
img(:,:,dim) = im;
end
rect(1) = rect(1)+shifts(1);
rect(1) = rect(1)-shifts(2);
rect(2) = rect(2)+shifts(3);
rect(2) = rect(2)-shifts(4);

img = imcrop(img,rect);
end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step6_v1_4(specs,show)
%% CreateFishStruct_step6
%  Version 1.4
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 2/12/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Use the rough mask/outline (R_body) of fish body from earlier to find better more accurate outline (body)
%  save in rounds(2)
%% Version History
%  1.1: using return_body_outlines_mod2; changed how fin bones are removed
%  and how tail position is set
%  1.2: processing rgb images
%  1.3: return_body_outlines_mod5
%  1.4: one fish had an incorrect body outline that included tail. Trying
%  to fix (using return_body_outlines_mod6)

disp('Step 6: Finding more accurate outline of fish bodies and save in rounds(2)...');

tic;

for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}
% retrieave images from data struct.
if(ndims(FISH(ff).rounds(2).bf1)==3)
% red channel is best because tail doesn't show up as much
bf1   = FISH(ff).rounds(2).bf1(:,:,1);
bf2   = FISH(ff).rounds(2).bf2(:,:,1);
%             bf1   = rgb2gray(FISH(ff).rounds(2).bf1);
%             bf2   = rgb2gray(FISH(ff).rounds(2).bf2);
else
bf1   = FISH(ff).rounds(2).bf1;
bf2   = FISH(ff).rounds(2).bf2;
end
rgbr1 = FISH(ff).rounds(2).rgb1(:,:,1);
rgbr2 = FISH(ff).rounds(2).rgb2(:,:,1);

R_body1 = FISH(ff).rounds(2).R_body1;
R_body2 = FISH(ff).rounds(2).R_body2;
R_body3 = FISH(ff).rounds(2).R_body3;
R_body4 = FISH(ff).rounds(2).R_body4;

% the function 'return_body_outlines' uses old rough outline
% R_body and the bf images of both sides 
% to find new more accurate boundaries of 
% the fish bodies
[body1,body2,body1_end,body2_end,edge1,edge2,edge_box1,edge_box2,...
shift1] = return_body_outlines_mod6(bf1,R_body1,bf2,R_body2,specs.minimum_area);
[body3,body4,body3_end,body4_end,c_edge1,c_edge2,c_edge_box1,...
c_edge_box2,cshift1] = return_body_outlines_mod6(rgbr1,R_body3,rgbr2,R_body4,specs.minimum_area);

% determine points along boundary for plotting
B1 = bwboundaries(body1);
B2 = bwboundaries(body2);
B3 = bwboundaries(body3);
B4 = bwboundaries(body4);

if show==1
% visually inspect the new improved outlines plotted on top of BF images
figure(1)
subplot(2,2,1)
imshow(bf1)
hold on
for k=1:length(B1)
b=B1{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end
title(['Fish number: ' num2str(ff) ' ,  Time: ' num2str(t) '  - Right side' ],'FontSize',10)
subplot(2,2,3)
imshow(bf2)
hold on
for k=1:length(B2)
b=B2{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end        
title(['Fish number: ' num2str(ff) '  ,  Time: ' num2str(t) '  - Left side' ],'FontSize',10)
subplot(2,2,2)
imshow(rgbr1)
hold on
for k=1:length(B3)
b=B3{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end
title(['Fish number: ' num2str(ff) ' ,  Time: ' num2str(t) '  - Right side' ],'FontSize',10)
subplot(2,2,4)
imshow(rgbr2)
hold on
for k=1:length(B4)
b=B4{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end        
title(['Fish number: ' num2str(ff) '  ,  Time: ' num2str(t) '  - Left side' ],'FontSize',10)
pause(0.1);
end

% Save improved body outline in rounds(2)
FISH(ff).rounds(2).body1 = body1;
FISH(ff).rounds(2).body2 = body2;
FISH(ff).rounds(2).body3 = body3;
FISH(ff).rounds(2).body4 = body4;
% save end points for body
FISH(ff).rounds(2).bodyEnd1_coor = body1_end;
FISH(ff).rounds(2).bodyEnd2_coor = body2_end;
FISH(ff).rounds(2).bodyEnd3_coor = body3_end;
FISH(ff).rounds(2).bodyEnd4_coor = body4_end;
% save tail edge images
FISH(ff).rounds(2).tail_edge1 = edge1;
FISH(ff).rounds(2).tail_edge2 = edge2;
FISH(ff).rounds(2).tail_edge3 = c_edge1;
FISH(ff).rounds(2).tail_edge4 = c_edge2;
FISH(ff).rounds(2).tail_box1 = edge_box1;
FISH(ff).rounds(2).tail_box2 = edge_box2;
FISH(ff).rounds(2).tail_box3 = c_edge_box1;
FISH(ff).rounds(2).tail_box4 = c_edge_box2;
FISH(ff).rounds(2).tail_LR_shift_bf = shift1;
FISH(ff).rounds(2).tail_LR_shift_rgb = cshift1;
end
BYTIME(t).fishes = FISH;
end
toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step8_v1_3(specs,show)
%% CreateFishStruct_step8
%  Version 1.3
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 9/25/19
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Find and save eye outlines and center of mass (landmark) in rounds(2)
%  1.2: process rgb images 
%  1.3: want to refine eye recognition

disp('Step 8: Finding eye outlines and eye center of mass landmark and save in rounds(2)...');

tic;

% put in approx coordinates of eye and approx diameter of eye
% approx_eye_coor = [75 180];
approx_eye_coor = specs.eye_coordinate;
approx_eye_diameter = specs.eye_diameter;
CM_eye_avg = [0 0];
avg_eye_diamter = 0;
nfish = 0;
switch_text = {'BF1','BF2','RGB1','RGB2'};
% rectangle for cropping out part of eye, used by function
% 'return_eye_outline' inside loop
% rect = [approx_eye_coor(1)-x_shift     approx_eye_coor(2)-y_shift     approx_eye_diameter*5/3     approx_eye_diameter/2]; % [xmax ymax width height]
% rect = [approx_eye_coor(1)-x_shift     approx_eye_coor(2)-y_shift     x_shift*5     y_shift*2]; % [xmax ymax width height]

for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}
LRLR = 0;

while LRLR<4
LRLR = LRLR+1;
% retrieave images from data struct.
switch LRLR
case 1
bf  = FISH(ff).rounds(2).bf1(:,:,1);
body = FISH(ff).rounds(2).body1;
itsf = specs.eyeBrightnessFactor; % makes dark eyes stand out
case 2
bf  = FISH(ff).rounds(2).bf2(:,:,1);
body = FISH(ff).rounds(2).body2;
itsf = specs.eyeBrightnessFactor; % makes dark eyes stand out
case 3
bf = FISH(ff).rounds(2).rgb1(:,:,1);
body = FISH(ff).rounds(2).body3;
itsf = 0.33*specs.eyeBrightnessFactor; % makes dark eyes stand out
case 4
bf = FISH(ff).rounds(2).rgb2(:,:,1);
body = FISH(ff).rounds(2).body4;
itsf = 0.33*specs.eyeBrightnessFactor; % makes dark eyes stand out
end

% find eye circle (hopefully)
crps = 31;
crpe = approx_eye_coor(1) + approx_eye_diameter;
row_range = approx_eye_coor(2)+(-150:150);
% get better itsf value
bfhead = bf(row_range,crps:crpe);
bfhead(body(row_range,crps:crpe)==0)=0;
itsf2 = 65535/double(median(nonzeros(bfhead)));
itsf = min(itsf,itsf2);
[cents,rads] = imfindcircles(itsf*bf(row_range,crps:crpe),round(0.5*approx_eye_diameter)+[-7 7],...
'ObjectPolarity','dark','Sensitivity',0.94,'EdgeThreshold',0.01);
if(isempty(cents))
real_eye_coor = approx_eye_coor;
real_eye_diameter = approx_eye_diameter;
elseif(length(cents(:,1))==1)
cents(:,2) = cents(:,2) + row_range(1)-1;
cents = cents + [crps-1 0];
real_eye_coor = round(cents);
real_eye_diameter = round(2*rads);
else
% peak circle with centroid closest to approx_eye_coor
cents(:,2) = cents(:,2) + row_range(1)-1;
cents = cents + [crps-1 0];
dists = (cents(:,1)-approx_eye_coor(1)).^2+(cents(:,2)-approx_eye_coor(2)).^2;
[~,pick] = min(dists);
real_eye_coor = round(cents(pick,:));
real_eye_diameter = round(2*rads(pick));
end
x_shift = round(real_eye_diameter/3);
y_shift = round(real_eye_diameter/4);
rect = [real_eye_coor(1)-x_shift     real_eye_coor(2)-y_shift     x_shift*5     y_shift*2]; % [xmax ymax width height]

% function 'return_eye_outline' takes the image
% and rectangle (based on approx. eye center and diameter) and return
% eye outline/mask
max_crop = double(max(max(bf(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3)))));
itsf3 = max(1,min(5e4/max_crop,itsf));
bw = return_eye_outline_mod4(itsf3*bf,rect);
% create a circle mask if imfindcircles has a result
if(~isempty(cents))
[row,col] = size(bw);
xdv = (1:col)-real_eye_coor(1);
ydv = (1:row)-real_eye_coor(2);
[xds,yds] = meshgrid(xdv.^2,ydv.^2);
dist_from_eye_sqr = xds+yds;
bw_circ1 = dist_from_eye_sqr <= (real_eye_diameter*0.48)^2;
bw(bw_circ1==1) = 1;
bw_circ2 = dist_from_eye_sqr <= (real_eye_diameter*0.55)^2;
bw(bw_circ2==0) = 0;
end

% if eye not detected, then fix it
eye_area = sum(bw(:));
max_area = (real_eye_diameter+20)^2*pi()/4;
if(eye_area < specs.minimum_eye_area)
% put in a circle based on approx eye location
[row,col] = size(bw);
if(eye_area > 0)
[yeye,xeye] = find(bw==1);
est_eye_coor = [round(mean(xeye)) round(mean(yeye))];
else
est_eye_coor = approx_eye_coor;
end
xdv = (1:col)-est_eye_coor(1);
ydv = (1:row)-est_eye_coor(2);
[xds,yds] = meshgrid(xdv.^2,ydv.^2);
dist_from_eye_sqr = xds+yds;
bw_circ1 = dist_from_eye_sqr <= (approx_eye_diameter*0.4)^2;
bw(bw_circ1==1) = 1;
warning('No eye detected for fish: %i at time: %i. (Approximate one was drawn)',ff,t);
elseif(eye_area > max_area)
warning('Eye is too large for fish: %i at time: %i (%f>%f)',ff,t,eye_area,max_area);
end

r = regionprops(bw,'Centroid');
if(length(r) > 1)
warning('Too many eyes in outline at fish %i, time %i. (Selected largest one)',ff,t);
r2 = regionprops(bw,'Area');
imax_area = 1;
for ir = 2:length(r2)
if(r2(ir).Area > r2(imax_area).Area)
imax_area = ir;
end
end
bwlab = bwlabel(bw);
bw = false(size(bw));
bw(bwlab==imax_area) = 1;
r = regionprops(bw,'Centroid');
end

CM_eye = [r.Centroid];
CM_eye_avg = CM_eye_avg + CM_eye;
avg_eye_diamter = avg_eye_diamter + real_eye_diameter;
nfish = nfish + 1;

if show == 1
% visually inspect that eye outline was determined correctly
B = bwboundaries(bw);
figure(1)
imshow(itsf*bf)
hold on
for k=1:length(B)
b=B{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end
plot(CM_eye(1),CM_eye(2),'g*','MarkerSize',15)
plot([rect(1) rect(1) rect(1)+rect(3) rect(1)+rect(3) rect(1)], [rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2) rect(2)],'b--');
if(~isempty(cents))
viscircles(real_eye_coor,0.5*real_eye_diameter,'LineStyle',':','LineWidth',1);
end
title(['Eye outline and center ' switch_text{LRLR} ' - Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11)
pause(0.1);
end

% save eye outline and CM coor in rounds(2)
switch LRLR
case 1
FISH(ff).rounds(2).eye1 = bw;
FISH(ff).rounds(2).eye1_coor = round(CM_eye);
case 2
FISH(ff).rounds(2).eye2 = bw;
FISH(ff).rounds(2).eye2_coor = round(CM_eye);
case 3
FISH(ff).rounds(2).eye3 = bw;
FISH(ff).rounds(2).eye3_coor = round(CM_eye);
case 4
FISH(ff).rounds(2).eye4 = bw;
FISH(ff).rounds(2).eye4_coor = round(CM_eye);

end

end
end
BYTIME(t).fishes = FISH;
end

CM_eye_avg = CM_eye_avg/nfish;
avg_eye_diamter = avg_eye_diamter/nfish;
fprintf('Average eye position: (%4.1f, %4.1f)\n',CM_eye_avg(1),CM_eye_avg(2));
fprintf('Average eye diamter: %4.1f\n',avg_eye_diamter);
toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step7_v2_2(name,show)
%% CreateFishStruct_step7
%  Version 2.2
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 11/15/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Find spine outline and end of spine landmark point
%% Version History
%  1.1: changed to return_spine_outline_mod
%  2.0: changed to return_spine_outline_mod3 (better spine detection and
%  also returns a polynomial fit of spine) Also uses eye position and body
%  end to get better spine function
%  2.1: processing rgb images
%  2.2: return_spine_outline_mod4

disp('Step 7: Finding spine outlines and end of spine landmark point and save in rounds(2)...');

tic;
if(~isdir(name))
mkdir(name);
end
if(~isdir([name '/spines']))
mkdir([name '/spines']);
end
for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}
% retrieave images from data struct.
bf1   = FISH(ff).rounds(2).bf1;
body1 = FISH(ff).rounds(2).body1;
gfp1 = FISH(ff).rounds(2).gfp1;
bf2   = FISH(ff).rounds(2).bf2;
body2 = FISH(ff).rounds(2).body2;
gfp2 = FISH(ff).rounds(2).gfp2;
rgb1 = FISH(ff).rounds(2).rgb1;
body3 = FISH(ff).rounds(2).body3;
rgb2 = FISH(ff).rounds(2).rgb2;
body4 = FISH(ff).rounds(2).body4;

eye1_coor = FISH(ff).rounds(2).eye1_coor;
eye2_coor = FISH(ff).rounds(2).eye2_coor;
body1_end = FISH(ff).rounds(2).bodyEnd1_coor;
body2_end = FISH(ff).rounds(2).bodyEnd2_coor;
eye3_coor = FISH(ff).rounds(2).eye3_coor;
eye4_coor = FISH(ff).rounds(2).eye4_coor;
body3_end = FISH(ff).rounds(2).bodyEnd3_coor;
body4_end = FISH(ff).rounds(2).bodyEnd4_coor;

% function 'return_spine_outline' takes red channel of bf
% and body outline/mask and return spine outline/mask
out_name1 = [name '/spines/' FISH(ff).fishID '_spine1'];
out_name2 = [name '/spines/' FISH(ff).fishID '_spine2'];
out_name3 = [name '/spines/' FISH(ff).fishID '_spine3'];
out_name4 = [name '/spines/' FISH(ff).fishID '_spine4'];
[spine1,spine_func1] = return_spine_outline_mod5(bf1(:,:,1),body1,body1_end,eye1_coor,gfp1,out_name1);
[spine2,spine_func2] = return_spine_outline_mod5(bf2(:,:,1),body2,body2_end,eye2_coor,gfp2,out_name2);
[spine3,cspine_func1] = return_spine_outline_mod5(rgb1(:,:,1),body3,body3_end,eye3_coor,gfp1,out_name3);
[spine4,cspine_func2] = return_spine_outline_mod5(rgb2(:,:,1),body4,body4_end,eye4_coor,gfp2,out_name4);

% get spine end points
spine1s = smoothBW_mod(spine1,20);
spine2s = smoothBW_mod(spine2,20);
spine3s = smoothBW_mod(spine3,20);
spine4s = smoothBW_mod(spine4,20);

spine1_end = getBodyEndPoint_v1_0(spine1s,20,70);
spine2_end = getBodyEndPoint_v1_0(spine2s,20,70);
spine3_end = getBodyEndPoint_v1_0(spine3s,20,70);
spine4_end = getBodyEndPoint_v1_0(spine4s,20,70);

if show==1
% plot spine outlines and
% end of spine landmarks to visually inspect that they look right
S1 = bwboundaries(spine1);
S2 = bwboundaries(spine2); 
S3 = bwboundaries(spine3);
S4 = bwboundaries(spine4); 
figure(1)
subplot(2,1,1)
imshowpair(bf1,bf2)
hold on
plot(body1_end(1),body1_end(2),'b*','LineWidth',2);
plot(body2_end(1),body2_end(2),'r*','LineWidth',2);
plot(spine1_end(1),spine1_end(2),'bo','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','y')
plot(spine2_end(1),spine2_end(2),'ro','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','y')
for k=1:length(S1)
b=S1{k};
plot(b(:,2),b(:,1),'b-','LineWidth',2)
end
for k=1:length(S2)
b=S2{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end
title(['Spine - Fish: ' num2str(ff) ',  Time: ' num2str(t) ' , Rigth side blue, Left side red'],'FontSize',11)
subplot(2,1,2)
imshowpair(rgb1,rgb2)
hold on
plot(body3_end(1),body3_end(2),'b*','LineWidth',2);
plot(body4_end(1),body4_end(2),'r*','LineWidth',2);
plot(spine3_end(1),spine3_end(2),'bo','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','y')
plot(spine4_end(1),spine4_end(2),'ro','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','y')
for k=1:length(S3)
b=S3{k};
plot(b(:,2),b(:,1),'b-','LineWidth',2)
end
for k=1:length(S4)
b=S4{k};
plot(b(:,2),b(:,1),'r-','LineWidth',2)
end
title(['Spine - Fish: ' num2str(ff) ',  Time: ' num2str(t) ' , Rigth side blue, Left side red'],'FontSize',11)

pause(0.1);
end

% save spine outlines in rounds(2)
FISH(ff).rounds(2).spine1 = spine1;
FISH(ff).rounds(2).spine2 = spine2;
FISH(ff).rounds(2).spine_details1 = spine_func1;
FISH(ff).rounds(2).spine_details2 = spine_func2;
FISH(ff).rounds(2).spine3 = spine3;
FISH(ff).rounds(2).spine4 = spine4;
FISH(ff).rounds(2).spine_details3 = cspine_func1;
FISH(ff).rounds(2).spine_details4 = cspine_func2;

% save end of spine landmark coor in rounds(2)
FISH(ff).rounds(2).spineEnd1_coor = spine1_end;
FISH(ff).rounds(2).spineEnd2_coor = spine2_end;
FISH(ff).rounds(2).spineEnd3_coor = spine3_end;
FISH(ff).rounds(2).spineEnd4_coor = spine4_end;

end
BYTIME(t).fishes = FISH;
end
toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step9_v2_1(specs,show)
%% CreateFishStruct_step9
%  Version 2.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 10/20/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Find boundary landmark points and save left-rigth transforms in rounds(2)
% transform images and save transformed images in rounds(3)
%
%% Version History
%  1.1: using body end as alignment point instead of spine end (this is
%  after the body ends have been pruned flat (CreateFishStruct_step6_v1_1)
%  changed ctp2trans and imtranfrom to fitgeotrans and imwarp
%  2.0: need separate x_coor for each side of the fish; otherwise the two
%  sides are different lengths after transformation
%  2.1: processing rgb images

disp('Step 9: Finding boundary landmark points... doing left-rigth side transforms... saving transformed images in rounds(3)...' );

tic;

% get width and heigth of images
ww = specs.maximum_width;
hh = specs.maximum_length;
objs(2) = struct('spine','input_points');
objs2(2) = struct('spine','input_points');
for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}

% retrive images from data struc
body1  = FISH(ff).rounds(2).body1;
body2  = FISH(ff).rounds(2).body2;
eye1  = FISH(ff).rounds(2).eye1;
eye2  = FISH(ff).rounds(2).eye2;
spine1  = FISH(ff).rounds(2).spine1;
spine2  = FISH(ff).rounds(2).spine2;
bf1  = FISH(ff).rounds(2).bf1;
bf2  = FISH(ff).rounds(2).bf2;
gfp1  = FISH(ff).rounds(2).gfp1;
gfp2  = FISH(ff).rounds(2).gfp2;
rfp1  = FISH(ff).rounds(2).rfp1;
rfp2  = FISH(ff).rounds(2).rfp2;

body3 = FISH(ff).rounds(2).body3;
body4 = FISH(ff).rounds(2).body4;
eye3  = FISH(ff).rounds(2).eye3;
eye4  = FISH(ff).rounds(2).eye4;
spine3  = FISH(ff).rounds(2).spine3;
spine4  = FISH(ff).rounds(2).spine4;
rgb1 = FISH(ff).rounds(2).rgb1;
rgb2 = FISH(ff).rounds(2).rgb2;

% retrive landmark coor from data struc
eye_coor1 = round(FISH(ff).rounds(2).eye1_coor);
eye_coor2 = round(FISH(ff).rounds(2).eye2_coor);
end_coor1 = FISH(ff).rounds(2).bodyEnd1_coor;
end_coor2 = FISH(ff).rounds(2).bodyEnd2_coor;

eye_coor3 = round(FISH(ff).rounds(2).eye3_coor);
eye_coor4 = round(FISH(ff).rounds(2).eye4_coor);
end_coor3 = FISH(ff).rounds(2).bodyEnd3_coor;
end_coor4 = FISH(ff).rounds(2).bodyEnd4_coor;

%         if(size(eye_coor1)~=size(eye_coor2))
%             disp('stop');
%         end
% make average of eye and spine landmark point the base points
eye_coor_base = round((eye_coor1+eye_coor2)./2);
end_coor_base = round((end_coor1+end_coor2)./2);
eye_coor_base2 = round((eye_coor3+eye_coor4)./2);
end_coor_base2 = round((end_coor3+end_coor4)./2);

% use spine to get reference points for transformation
objs(1).spine = spine1;
objs(2).spine = spine2;
objs(1).input_points = [eye_coor1;end_coor1];
objs(2).input_points = [eye_coor2;end_coor2];

objs2(1).spine = spine3;
objs2(2).spine = spine4;
objs2(1).input_points = [eye_coor3;end_coor3];
objs2(2).input_points = [eye_coor4;end_coor4];

[x_coor_spine,y_coor_spine,x_coor_spine_base,y_coor_spine_base,spine_slope] = return_spine_coor_T_mod2(objs,6);
[x_coor_spine2,y_coor_spine2,x_coor_spine_base2,y_coor_spine_base2,spine_slope2] = return_spine_coor_T_mod2(objs2,6);

[TFORM1,TFORM2,input_points1,input_points2,base_points] = getTransformation(...
eye_coor1,eye_coor2,eye_coor_base,end_coor1,end_coor2,end_coor_base,...
x_coor_spine,y_coor_spine,x_coor_spine_base,y_coor_spine_base,spine_slope);
[TFORM3,TFORM4,input_points3,input_points4,base_points2] = getTransformation(...
eye_coor3,eye_coor4,eye_coor_base2,end_coor3,end_coor4,end_coor_base2,...
x_coor_spine2,y_coor_spine2,x_coor_spine_base2,y_coor_spine_base2,spine_slope2);
%         % concatenate x- and y- coor of input and base points
%         x_coor1 = [x_coor_spine(1,:) eye_coor1(1) end_coor1(1)];
%         x_coor2 = [x_coor_spine(2,:) eye_coor2(1) end_coor2(1)];
%         x_coor_base = [x_coor_spine_base eye_coor_base(1) end_coor_base(1)];
%         y_coor1 = [y_coor_spine(1,:) eye_coor1(2) end_coor1(2)];
%         y_coor2 = [y_coor_spine(2,:) eye_coor2(2) end_coor2(2)];
%         y_coor_base = [y_coor_spine_base eye_coor_base(2) end_coor_base(2)];
%         
%         
%         % other points to help with transformation - these points define
%         % alignment in the y direction - all the points are spaced the same
%         % amount so that mostly translation and not stretching occurs in
%         % the y direction.
%         up_spacing = 150;
%         down_spacing = 200;
%         y_up1 = y_coor1-up_spacing;
%         y_down1 = y_coor1+down_spacing;
%         y_up2 = y_coor2-up_spacing;
%         y_down2 = y_coor2+down_spacing;
%         y_up_base = y_coor_base-up_spacing;
%         y_down_base = y_coor_base+down_spacing;
%         
%         
%         
%         % initialize input point and base point vectors
%         input_points1 = zeros(length(x_coor1)*3,2);
%         input_points2 = zeros(length(x_coor2)*3,2);
%         base_points   = zeros(length(x_coor_base)*3,2);
%         
%         % fill in values
%         input_points1(:,1) = [x_coor1 x_coor1 x_coor1]'; 
%         input_points2(:,1) = [x_coor2 x_coor2 x_coor2]'; 
%         base_points(:,1) = [x_coor_base x_coor_base x_coor_base]';
%         input_points1(:,2) = [y_up1 y_coor1 y_down1]'; 
%         input_points2(:,2) = [y_up2 y_coor2 y_down2]'; 
%         base_points(:,2) = [y_up_base y_coor_base y_down_base]';
%         
%         % calculate image tranformation functions
%         TFORM1 = fitgeotrans(input_points1, base_points,'polynomial',2);
%         TFORM2 = fitgeotrans(input_points2, base_points,'polynomial',2);
% get imref2d object
RA = imref2d(size(body1),[0.5 hh+0.5],[0.5 ww+0.5]);
bf_fill = round(mean([bf1(bf1>mean(bf1(:))); bf2(bf2>mean(bf2(:)))]));
% apply tranforms to all images
bodyT1 = imwarp(body1,TFORM1,'OutputView',RA,'FillValues',0);
bodyT2 = imwarp(body2,TFORM2,'OutputView',RA,'FillValues',0);
spineT1 = imwarp(spine1,TFORM1,'OutputView',RA,'FillValues',0);
spineT2 = imwarp(spine2,TFORM2,'OutputView',RA,'FillValues',0);
bfT1 = imwarp(bf1,TFORM1,'OutputView',RA,'FillValues',bf_fill);
bfT2 = imwarp(bf2,TFORM2,'OutputView',RA,'FillValues',bf_fill);
eyeT1 = imwarp(eye1,TFORM1,'OutputView',RA,'FillValues',0);
gfpT1 = imwarp(gfp1,TFORM1,'OutputView',RA,'FillValues',0);
rfpT1 = imwarp(rfp1,TFORM1,'OutputView',RA,'FillValues',0);
eyeT2 = imwarp(eye2,TFORM2,'OutputView',RA,'FillValues',0);
gfpT2 = imwarp(gfp2,TFORM2,'OutputView',RA,'FillValues',0);
rfpT2 = imwarp(rfp2,TFORM2,'OutputView',RA,'FillValues',0);

% warp cbody, eye and all three spectra or rgb image
bodyT3 = imwarp(body3,TFORM3,'OutputView',RA,'FillValues',0);
bodyT4 = imwarp(body4,TFORM4,'OutputView',RA,'FillValues',0);
spineT3 = imwarp(spine3,TFORM3,'OutputView',RA,'FillValues',0);
spineT4 = imwarp(spine4,TFORM4,'OutputView',RA,'FillValues',0);
eyeT3 = imwarp(eye3,TFORM3,'OutputView',RA,'FillValues',0);
eyeT4 = imwarp(eye4,TFORM4,'OutputView',RA,'FillValues',0);
rgbT1 = rgb1;
rgbT2 = rgb2;
for dim = 1:size(rgb1,3)
temp = rgb1(:,:,dim);
temp_fill = round(mean(temp(temp>mean(temp(:)))));
tempT = imwarp(temp,TFORM3,'OutputView',RA,'FillValues',temp_fill);
rgbT1(:,:,dim) = tempT;
end
for dim = 1:size(rgb2,3)
temp = rgb2(:,:,dim);
temp_fill = round(mean(temp(temp>mean(temp(:)))));
tempT = imwarp(temp,TFORM4,'OutputView',RA,'FillValues',temp_fill);
rgbT2(:,:,dim) = tempT;
end

% save left-right side transformed images in rounds(3)
FISH(ff).rounds(3).bf1   = uint16(bfT1);
FISH(ff).rounds(3).bf2   = uint16(bfT2);
FISH(ff).rounds(3).gfp1  = uint16(gfpT1);
FISH(ff).rounds(3).gfp2  = uint16(gfpT2);
FISH(ff).rounds(3).rfp1  = uint16(rfpT1);
FISH(ff).rounds(3).rfp2  = uint16(rfpT2);
FISH(ff).rounds(3).eye1  = logical(eyeT1);
FISH(ff).rounds(3).eye2  = logical(eyeT2);
FISH(ff).rounds(3).body1 = logical(bodyT1);
FISH(ff).rounds(3).body2 = logical(bodyT2);
FISH(ff).rounds(3).spine1 = logical(spineT1);
FISH(ff).rounds(3).spine2 = logical(spineT2);

FISH(ff).rounds(3).rgb1 = uint16(rgbT1);
FISH(ff).rounds(3).rgb2 = uint16(rgbT2);
FISH(ff).rounds(3).eye3 = uint16(eyeT3);
FISH(ff).rounds(3).eye4 = uint16(eyeT4);
FISH(ff).rounds(3).body3 = uint16(bodyT3);
FISH(ff).rounds(3).body4 = uint16(bodyT4);
FISH(ff).rounds(3).spine3 = uint16(spineT3);
FISH(ff).rounds(3).spine4 = uint16(spineT4);

% save transform functions        
FISH(ff).transforms.left_rightTF1 = TFORM1;
FISH(ff).transforms.left_rightTF1_distortion = sum(bodyT1(:))/sum(body1(:));
FISH(ff).transforms.left_rightTF2 = TFORM2;
FISH(ff).transforms.left_rightTF2_distortion = sum(bodyT2(:))/sum(body2(:));
FISH(ff).transforms.left_rightTF3 = TFORM3;
FISH(ff).transforms.left_rightTF3_distortion = sum(bodyT3(:))/sum(body3(:));
FISH(ff).transforms.left_rightTF4 = TFORM4;
FISH(ff).transforms.left_rightTF4_distortion = sum(bodyT4(:))/sum(body4(:));


% get new end positions
% guess is the new position if only translation occured
guess_coor1 = 2*end_coor1-getTransformPointsAsArray(TFORM1,end_coor1);
guess_coor2 = 2*end_coor2-getTransformPointsAsArray(TFORM2,end_coor2);
guess_coor3 = 2*end_coor3-getTransformPointsAsArray(TFORM3,end_coor3);
guess_coor4 = 2*end_coor4-getTransformPointsAsArray(TFORM4,end_coor4);
% will use optimization to find accurate solution
opt_func1 = @(x) sumsqr(getTransformPointsAsArray(TFORM1,x)-end_coor1);
opt_func2 = @(x) sumsqr(getTransformPointsAsArray(TFORM2,x)-end_coor2);
opt_func3 = @(x) sumsqr(getTransformPointsAsArray(TFORM3,x)-end_coor3);
opt_func4 = @(x) sumsqr(getTransformPointsAsArray(TFORM4,x)-end_coor4);
new_coor1 = round(fminsearch(opt_func1,guess_coor1));
new_coor2 = round(fminsearch(opt_func2,guess_coor2));
new_coor3 = round(fminsearch(opt_func3,guess_coor3));
new_coor4 = round(fminsearch(opt_func4,guess_coor4));

FISH(ff).rounds(3).bodyEnd1_coor = new_coor1;
FISH(ff).rounds(3).bodyEnd2_coor = new_coor2;
FISH(ff).rounds(3).bodyEnd3_coor = new_coor3;
FISH(ff).rounds(3).bodyEnd4_coor = new_coor4;

% get new eye positons
r1 = regionprops(FISH(ff).rounds(3).eye1,'Centroid');
FISH(ff).rounds(3).eye1_coor = r1.Centroid;
r2 = regionprops(FISH(ff).rounds(3).eye2,'Centroid');
FISH(ff).rounds(3).eye2_coor = r2.Centroid;
r3 = regionprops(FISH(ff).rounds(3).eye3,'Centroid');
FISH(ff).rounds(3).eye3_coor = r3.Centroid;
r4 = regionprops(FISH(ff).rounds(3).eye4,'Centroid');
FISH(ff).rounds(3).eye4_coor = r4.Centroid;

if(length(r1)~=1 || length(r2)~=1 || length(r3)~=1 || length(r4)~=1)
error('Bad eye outline for fish %i, time %i',ff,t);
end

% visually inspect if transform
% make left and right side fish bodies overlay nicely 
if show==1
bf1_show = bf1;
bf1_show(body1==0)=0;
bf2_show = bf2;
bf2_show(body2==0)=0;
figure(1)
subplot(2,2,1); 
imshowpair(bf1_show,bf2_show);
hold on;
A1 = bwboundaries(body1-eye1);
for k=1:length(A1)
b=A1{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
A2 = bwboundaries(body2-eye2);
for k=1:length(A2)
b=A2{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
plot(input_points1(:,1),input_points1(:,2),'+r','MarkerSize',15,'LineWidth',1)
plot(input_points2(:,1),input_points2(:,2),'xb','MarkerSize',15,'LineWidth',1)
plot(base_points(:,1),base_points(:,2),'og','MarkerSize',15,'LineWidth',1)
%             legend('Right side input points','Left side input points','Base points')
title(['Originals - Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11);
hold off;
bfT1_show = bfT1;
bfT1_show(bodyT1==0)=0;
bfT2_show = bfT2;
bfT2_show(bodyT2==0)=0;
subplot(2,2,3); 
imshowpair(bfT1_show,bfT2_show);
hold on;
B1 = bwboundaries(bodyT1-eyeT1);
for k=1:length(B1)
b=B1{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
B2 = bwboundaries(bodyT2-eyeT2);
for k=1:length(B2)
b=B2{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
B3 = bwboundaries(spineT1);
for k=1:length(B3)
b=B3{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
B4 = bwboundaries(spineT2);
for k=1:length(B4)
b=B4{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
plot([new_coor1(1) new_coor2(1)],[new_coor1(2) new_coor2(2)],'og','MarkerSize',10,'LineWidth',1);
plot([r1.Centroid(1) r2.Centroid(1)],[r1.Centroid(2) r2.Centroid(2)],'og','MarkerSize',10,'LineWidth',1);
title(['Transformed - Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11);
hold off;
subplot(2,2,2);
imshowpair(gfp1,gfp2);
hold on;
for k=1:length(A1)
b=A1{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
for k=1:length(A2)
b=A2{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
title(['GFP Original -  Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11);
subplot(2,2,4);
imshowpair(gfpT1,gfpT2);
hold on;
for k=1:length(B1)
b=B1{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
for k=1:length(B2)
b=B2{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
title(['GFP Transformed -  Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11);

figure(2)
subplot(2,1,1); 
imshowpair(rgb1,rgb2);
hold on;
A1 = bwboundaries(body3-eye3);
for k=1:length(A1)
b=A1{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
A2 = bwboundaries(body4-eye4);
for k=1:length(A2)
b=A2{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
plot(input_points3(:,1),input_points3(:,2),'+r','MarkerSize',15,'LineWidth',1)
plot(input_points4(:,1),input_points4(:,2),'xb','MarkerSize',15,'LineWidth',1)
plot(base_points2(:,1),base_points2(:,2),'og','MarkerSize',15,'LineWidth',1)
%             legend('Right side input points','Left side input points','Base points')
title(['Originals - Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11);
hold off;
subplot(2,1,2); 
imshowpair(rgbT1,rgbT2);
hold on;
B1 = bwboundaries(bodyT3-eyeT3);
for k=1:length(B1)
b=B1{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
B2 = bwboundaries(bodyT4-eyeT4);
for k=1:length(B2)
b=B2{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
B3 = bwboundaries(spineT3);
for k=1:length(B3)
b=B3{k};
plot(b(:,2),b(:,1),'-g','LineWidth',1);
end
B4 = bwboundaries(spineT4);
for k=1:length(B4)
b=B4{k};
plot(b(:,2),b(:,1),'-m','LineWidth',1);
end
plot([new_coor3(1) new_coor4(1)],[new_coor3(2) new_coor4(2)],'og','MarkerSize',10,'LineWidth',1);
plot([r3.Centroid(1) r4.Centroid(1)],[r3.Centroid(2) r4.Centroid(2)],'og','MarkerSize',10,'LineWidth',1);
title(['Transformed - Fish: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',11);
hold off;
pause(0.1);
end

end
BYTIME(t).fishes = FISH;
end

toc; 
disp(' ')
end

function [TFORM1,TFORM2,input_points1,input_points2,base_points] = getTransformation(eye_coor1,eye_coor2,eye_coor_base,end_coor1,end_coor2,end_coor_base,...
x_coor_spine,y_coor_spine,x_coor_spine_base,y_coor_spine_base,spine_slope)
% use anchor points to get transformation
% concatenate x- and y- coor of input and base points
x_coor1 = [x_coor_spine(1,:) eye_coor1(1) end_coor1(1)];
x_coor2 = [x_coor_spine(2,:) eye_coor2(1) end_coor2(1)];
x_coor_base = [x_coor_spine_base eye_coor_base(1) end_coor_base(1)];
y_coor1 = [y_coor_spine(1,:) eye_coor1(2) end_coor1(2)];
y_coor2 = [y_coor_spine(2,:) eye_coor2(2) end_coor2(2)];
y_coor_base = [y_coor_spine_base eye_coor_base(2) end_coor_base(2)];


% other points to help with transformation - these points define
% alignment in the y direction - all the points are spaced the same
% amount so that mostly translation and not stretching occurs in
% the y direction.
up_spacing = 150;
down_spacing = 200;
y_up1 = y_coor1-up_spacing;
y_down1 = y_coor1+down_spacing;
y_up2 = y_coor2-up_spacing;
y_down2 = y_coor2+down_spacing;
y_up_base = y_coor_base-up_spacing;
y_down_base = y_coor_base+down_spacing;



% initialize input point and base point vectors
input_points1 = zeros(length(x_coor1)*3,2);
input_points2 = zeros(length(x_coor2)*3,2);
base_points   = zeros(length(x_coor_base)*3,2);

% fill in values
input_points1(:,1) = [x_coor1 x_coor1 x_coor1]';
input_points2(:,1) = [x_coor2 x_coor2 x_coor2]';
base_points(:,1) = [x_coor_base x_coor_base x_coor_base]';
input_points1(:,2) = [y_up1 y_coor1 y_down1]';
input_points2(:,2) = [y_up2 y_coor2 y_down2]';
base_points(:,2) = [y_up_base y_coor_base y_down_base]';

% calculate image tranformation functions
TFORM1 = fitgeotrans(input_points1, base_points,'polynomial',2);
TFORM2 = fitgeotrans(input_points2, base_points,'polynomial',2);
end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step10_v1_3()
%% CreateFishStruct_step10
%  Version 1.2
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 10/20/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Now that left and rigth side images are aligned
% Merge bf, spine, body and eye images and save in rounds(3)
%% Version History:
%  1.1: will calculate and alternate body image using AND operation to make
% merged body rather than OR operation - this means the merged body will be
% smaller than both left and right sides
%  1.2: processing rgb files
%  1.3: adding a merge score

disp('Step 10: Making merged versions of left/right side images and saving in rounds(3)...' );

tic;

for t = 1:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}
% retrive images from data struc
body1  = FISH(ff).rounds(3).body1;
body2  = FISH(ff).rounds(3).body2;
eye1  = FISH(ff).rounds(3).eye1;
eye2  = FISH(ff).rounds(3).eye2;
spine1  = FISH(ff).rounds(3).spine1;
spine2  = FISH(ff).rounds(3).spine2;
bf1  = FISH(ff).rounds(3).bf1;
bf2  = FISH(ff).rounds(3).bf2;
body_end1 = FISH(ff).rounds(3).bodyEnd1_coor;
body_end2 = FISH(ff).rounds(3).bodyEnd2_coor;

body3  = FISH(ff).rounds(3).body3;
body4  = FISH(ff).rounds(3).body4;
eye3  = FISH(ff).rounds(3).eye3;
eye4  = FISH(ff).rounds(3).eye4;
spine3  = FISH(ff).rounds(3).spine3;
spine4  = FISH(ff).rounds(3).spine4;
rgb1  = FISH(ff).rounds(3).rgb1;
rgb2  = FISH(ff).rounds(3).rgb2;
body_end3 = FISH(ff).rounds(3).bodyEnd3_coor;
body_end4 = FISH(ff).rounds(3).bodyEnd4_coor;

% Merge bf1 and bf2 (pick darkest pixel from each image)
bf = MergeBF(bf1,bf2);
% Merge rgb images
rgb = rgb1;
for dim = 1:size(rgb1,3)
temp1 = rgb1(:,:,dim);
temp2 = rgb2(:,:,dim);
temp_merge = MergeBF(temp1,temp2);
rgb(:,:,dim) = temp_merge;
end

% Merge body1, body2 and spine1, spine2  and eye1, eye2 (pick all white pixels from each image)
union_body = body1;
union_body(body2==1)=1;
inter_body = ones(size(body1));
inter_body(body1==0)=0;
inter_body(body2==0)=0;
spine = spine1;
spine(spine2==1)=1;        
eye = eye1;
eye(eye2==1)=1;        
body_end = round(mean([body_end1; body_end2]));
r = regionprops(logical(eye),'Centroid');
eye_coor = round(r.Centroid);

union_body2 = body3;
union_body2(body4==1)=1;
inter_body2 = ones(size(body3));
inter_body2(body3==0)=0;
inter_body2(body4==0)=0;
cspine = spine3;
cspine(spine4==1)=1;        
ceye = eye3;
ceye(eye4==1)=1;        
cbody_end = round(mean([body_end3; body_end4]));
r2 = regionprops(logical(ceye),'Centroid');
ceye_coor = round(r2.Centroid);

% save in rounds(3)
FISH(ff).rounds(3).bf   = uint16(bf);
FISH(ff).rounds(3).eye  = logical(eye);
FISH(ff).rounds(3).body = logical(union_body);
FISH(ff).rounds(3).body_intersect = logical(inter_body);
FISH(ff).transforms.LR_merge_score = sum(inter_body(:))/sum(union_body(:));
FISH(ff).rounds(3).spine = logical(spine);
FISH(ff).rounds(3).bodyEnd_coor = body_end;
FISH(ff).rounds(3).eye_coor = eye_coor;

FISH(ff).rounds(3).rgb   = uint16(rgb);
FISH(ff).rounds(3).ceye  = logical(ceye);
FISH(ff).rounds(3).cbody = logical(union_body2);
FISH(ff).rounds(3).cbody_intersect = logical(inter_body2);
FISH(ff).transforms.LR_merge_score2 = sum(inter_body2(:))/sum(union_body2(:));
FISH(ff).rounds(3).cspine = logical(cspine);
FISH(ff).rounds(3).cbodyEnd_coor = cbody_end;
FISH(ff).rounds(3).ceye_coor = ceye_coor;
end     
BYTIME(t).fishes = FISH;
end
toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_identify_v2_1(lookup_file,show)
%% CreateFishStruct_identify_v2_1
%  Version 2.1
%  Author: Adeyinka Lesi
%  Date: 9/30/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Compare fish at different time points
%  figure out which fish is which from key points
%  1.1: added looking at gfp for extra reliability
%  2.0: implementing dataUpTo, tracking fish that do not survive all time
%  points
%  2.1: added option of a lookup file containing a verified configuration
%  also added a new term to the gfp error to penalize big shifts

disp('Step 10.5: Compare and identify fish');
tic;

% area in which to search for tumor
tum_rect = [200 150 800 300];
% other tunable variables
min_signal = 200;
fit_max = 3;


% open lookup file if one is giving
if(~isempty(lookup_file))
lookups = dlmread(lookup_file);
[nfish_lkup,tmax_lkup] = size(lookups);
% check viability of lookup file
% 1) right number of fish at each time point
nfish_t = max(lookups);
right_nfish = true;
for t = 1:tmax_lkup
if(length(BYTIME(end).dataAt{t})~=nfish_t(t))
right_nfish = false;
tmax_lkup = t;
break;
end
end
% 2) no repeats in any column (at any given time)
no_repeats = true;
for t = 1:tmax_lkup
if(length(unique(nonzeros(lookups(:,t))))~=nfish_t(t))
no_repeats = false;
end
end
% 3) all fish are accounted for
right_total = (nfish_lkup == length(BYTIME(1).fishes));

viable_lookup_file = all([right_nfish no_repeats right_total]);
else
tmax_lkup = 1;
viable_lookup_file = false;
end


% still need to create config for backwards compatibility
% fish (f,t+1) is the same as (config{t}(f),t)
config = cell(1,length(BYTIME)-1);
if(viable_lookup_file)
for t = 1:tmax_lkup-1
% each col of lookups contains the full set of fish numbers and
% possibly some zeros; the row then identifies which fid the fish
% belongs to
[~,fid_temp] = sort(lookups(:,t+1)); % fids in order of fish_number
fid_list = fid_temp(end-nfish_t(t+1)+1:end); % remove indices of zeros
% now can use fid to get fish numbers at time t
config{t} = lookups(fid_list,t);
end
end

if(viable_lookup_file && tmax_lkup >= length(BYTIME))
disp('Viable lookup table provided for all time points');
else    
% matrix of identities
% defacto order data
len = zeros(length(BYTIME(1).fishes),length(BYTIME));
wid = len;
eye1 = len;
eye2 = len;

if(viable_lookup_file)
fprintf('Viable lookup table provided for t=1 through t=%i\n',tmax_lkup);
else
disp('No viable lookup table provided');
end
fprintf('Attempting to assign identities for t=%i through t=%i\n',tmax_lkup+1,length(BYTIME));

for t = tmax_lkup:length(BYTIME)
FISH = BYTIME(t).fishes;
for ff=BYTIME(end).dataAt{t}
body = FISH(ff).rounds(3).body;
eye1c = FISH(ff).rounds(2).eye1_coor;
eye2c = FISH(ff).rounds(2).eye2_coor;
[y,x] = find(body,10);
nose = [mean(x) mean(y)];
[y,x] = find(body,100,'last');
tail = [mean(x) mean(y)];
len(ff,t) = tail(1)-nose(1);

eyec = 0.5*(eye1c+eye2c);
eye1(ff,t) = eyec(1);
eye2(ff,t) = eyec(2);

tbody = body';
[x,y] = find(tbody,100);
dors = [mean(x) mean(y)];
[x,y] = find(tbody,100,'last');
vent = [mean(x) mean(y)];
wid(ff,t) = vent(2)-dors(2);
end
end

% calculate overlap
overs = cell(1,length(BYTIME)-1);
lens = cell(1,length(BYTIME)-1);
wids = cell(1,length(BYTIME)-1);
eyes = cell(1,length(BYTIME)-1);
tots = cell(1,length(BYTIME)-1);
areas = cell(1,length(BYTIME)-1);
peris = cell(1,length(BYTIME)-1);
gfpfit1 = cell(1,length(BYTIME)-1);
gfpfit2 = cell(1,length(BYTIME)-1);
% tots{t}(f2,f) = [overs{t}(f2,f) lens{t}(f2,f) wids{t}(f2,f) eyes{t}(f2,f)]*weights;
weights = [1 1 2 4 0 0 4 4]';
weights = weights/sum(weights);

% to compare gfp images, need some preprocessing
for t = tmax_lkup:length(BYTIME)
for dai = 1:length(BYTIME(end).dataAt{t})
f = BYTIME(end).dataAt{t}(dai);
gfp1a = BYTIME(t).fishes(f).rounds(3).gfp1;
gfp2a = BYTIME(t).fishes(f).rounds(3).gfp2;
%         rfp1a = BYTIME(t).fishes(f).rounds(3).rfp1;
%         rfp2a = BYTIME(t).fishes(f).rounds(3).rfp2;
bod1a = BYTIME(t).fishes(f).rounds(3).body1;
bod2a = BYTIME(t).fishes(f).rounds(3).body2;
spn1a = BYTIME(t).fishes(f).rounds(3).spine1;
spn2a = BYTIME(t).fishes(f).rounds(3).spine2;
% adaptive threshold catches more stuff
glow1a = imbinarize(gfp1a,'adaptive');
glow2a = imbinarize(gfp2a,'adaptive');
% shape body outline to prune edges
[~,sc1] = ind2sub(size(spn1a),find(spn1a,1,'last'));
[~,sc2] = ind2sub(size(spn2a),find(spn2a,1,'last'));
bod1a2 = bod1a;
bod2a2 = bod2a;
bod1a2(:,sc1-100:end) = 0;
bod2a2(:,sc2-100:end) = 0;
bod1a2 = bwmorph(bod1a2,'erode',20);
bod2a2 = bwmorph(bod2a2,'erode',20);
% remove some noise by cropping the image
glow1a(bod1a2==0)=0;
glow2a(bod2a2==0)=0;
glow1a = bwmorph(glow1a,'erode',2);
glow2a = bwmorph(glow2a,'erode',2);
glow1a = bwmorph(glow1a,'dilate',2);
glow2a = bwmorph(glow2a,'dilate',2);

% draw outline
%         B1 = bwboundaries(bod1a);
%         B2 = bwboundaries(bod2a);
%         b_thickness = 5;
%         [img_ylen,img_xlen] = size(glow1a);
%         for c = 1:length(B1)
%             B = B1{c};
%             for bi = 1:size(B,1)
%                 % draw a box around each boundary point
%                 yhigh = min(img_ylen,round(B(bi,1) + 0.5*b_thickness));
%                 ylow = max(1,round(B(bi,1) - 0.5*b_thickness));
%                 xlow = max(1,round(B(bi,2) - 0.5*b_thickness));
%                 xhigh = min(img_xlen,round(B(bi,2) + 0.5*b_thickness));
%                 glow1a(ylow:yhigh,xlow:xhigh) = 1;
%             end
%         end
%         for c = 1:length(B2)
%             B = B2{c};
%             for bi = 1:size(B,1)
%                 % draw a box around each boundary point
%                 yhigh = min(img_ylen,round(B(bi,1) + 0.5*b_thickness));
%                 ylow = max(1,round(B(bi,1) - 0.5*b_thickness));
%                 xlow = max(1,round(B(bi,2) - 0.5*b_thickness));
%                 xhigh = min(img_xlen,round(B(bi,2) + 0.5*b_thickness));
%                 glow2a(ylow:yhigh,xlow:xhigh) = 1;
%             end
%         end

% to crop or not to crop
cglow1a = imcrop(glow1a,tum_rect);
cglow2a = imcrop(glow2a,tum_rect);
% save results
BYTIME(t).fishes(f).rounds(3).bw_gfp1 = glow1a;
BYTIME(t).fishes(f).rounds(3).bw_gfp2 = glow2a;
BYTIME(t).fishes(f).rounds(3).bw_gfp1_crop = cglow1a;
BYTIME(t).fishes(f).rounds(3).bw_gfp2_crop = cglow2a;
%         subplot(2,1,1);
%         imshow(double(bod1a)+double(bod1a2)+double(glow1a),[]);
%         title(sprintf('Side 1: t = %i, f = %i',t,f));
%         subplot(2,1,2);
%         imshow(double(bod2a)+double(bod2a2)+double(glow2a),[]);
%         title(sprintf('Side 2: t = %i, f = %i',t,f));
%         pause(0.1);
end
end

for t = tmax_lkup:length(BYTIME)-1
FISH = BYTIME(t).fishes;
FISH2 = BYTIME(t+1).fishes;
% array to measure interfish overlap compared to previous time
overs{t} = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
gfpfit1{t} = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
gfpfit2{t} = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
row_shift1 = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
row_shift2 = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
col_shift1 = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
col_shift2 = zeros(length(BYTIME(end).dataAt{t+1}),...
length(BYTIME(end).dataAt{t}));
for dai = 1:length(BYTIME(end).dataAt{t})
f = BYTIME(end).dataAt{t}(dai);
glow1a = FISH(f).rounds(3).bw_gfp1;
glow2a = FISH(f).rounds(3).bw_gfp2;
%         cglow1a = FISH(f).rounds(3).bw_gfp1_crop;
%         cglow2a = FISH(f).rounds(3).bw_gfp2_crop;
%         sum_glow1a = sum(cglow1a(:));
%         sum_glow2a = sum(cglow2a(:));
for dai2 = 1:length(BYTIME(end).dataAt{t+1})
f2 = BYTIME(end).dataAt{t+1}(dai2);
bod1 = FISH(f).rounds(3).body;
bod2 = FISH2(f2).rounds(3).body;
area1 = sum(bod1(:));
area2 = sum(bod2(:));
peri1 = sum(sum(bwperim(bod1)));
peri2 = sum(sum(bwperim(bod2)));
union = logical(bod1+bod2);
overlap = bod1;
overlap(bod2==0)=0;
nunn = sum(union(:));
novl = sum(overlap(:));
overs{t}(dai2,dai) = abs(log(novl/nunn));
areas{t}(dai2,dai) = abs(log(area1/area2));
peris{t}(dai2,dai) = abs(log(peri1/peri2));
lens{t}(dai2,dai) = abs(log(len(f,t)/len(f2,t+1)));
wids{t}(dai2,dai) = abs(log(wid(f,t)/wid(f2,t+1)));
eyes{t}(dai2,dai) = sqrt((eye1(f,t)-eye1(f2,t+1))^2+(eye2(f,t)-eye2(f2,t+1))^2)/len(f2,t+1);

cglow1b = FISH2(f2).rounds(3).bw_gfp1_crop;
cglow2b = FISH2(f2).rounds(3).bw_gfp2_crop;

% don't use gfp if their isn't a lot of signal
if(sum(cglow1b(:))>min_signal)
% align images (translates images over by a few pixels to
% maximize overlap)
%                 [sglow1a,dr1,dc1] = alignBW_v1_0(cglow1a,cglow1b);
[sglow1a,dr1,dc1] = alignBW_v2_0(glow1a,cglow1b,...
tum_rect(2),tum_rect(1),10,10);
row_shift1(dai2,dai) = dr1;
col_shift1(dai2,dai) = dc1;
% get overlap
ovlp1 = cglow1b;
ovlp1(sglow1a==0)=0;
unon1 = logical(sglow1a+cglow1b);
% calculate fit (will set a minimum for sum values to avoid
% errors)
sunn1 = max(sum(unon1(:)),1);
sovp1 = max(sum(ovlp1(:)),0.1); % 0.1 so that blank images don't have high fit
fit1 = abs(log(sovp1/sunn1))+(dr1^2+dc1^2)/800;
fit1 = min(fit1,fit_max);
%                 if(f==conf1(f2))
%                     figure(2);
%                     subplot(2,1,1);
%                     imshowpair(sglow1a,cglow1b);
%                     title(sprintf('(t%i,f%i) -> (t%i,f%i), fit: %3f',...
%                         t,f,t+1,f2,fit1));
%                 end
else
fit1 = 0;
end

if(sum(cglow2b(:))>min_signal)
% align images (translates images over by a few pixels to
% maximize overlap)
%                 [sglow2a,dr2,dc2] = alignBW_v1_0(cglow2a,cglow2b);
[sglow2a,dr2,dc2] = alignBW_v2_0(glow2a,cglow2b,...
tum_rect(2),tum_rect(1),10,10);
row_shift2(dai2,dai) = dr2;
col_shift2(dai2,dai) = dc2;
% get overlap
ovlp2 = cglow2b;
ovlp2(sglow2a==0)=0;
unon2 = logical(sglow2a+cglow2b);
% calculate fit (will set a minimum for sum values to avoid
% errors)
sunn2 = max(sum(unon2(:)),1);
sovp2 = max(sum(ovlp2(:)),0.1);
fit2 = abs(log(sovp2/sunn2))+(dr2^2+dc2^2)/800;
fit2 = min(fit2,fit_max);
%                 if(f==conf1(f2))
%                     figure(2);
%                     subplot(2,1,2);
%                     imshowpair(sglow2a,cglow2b);
%                     title(sprintf('(t%i,f%i) -> (t%i,f%i), fit: %3f',...
%                         t,f,t+1,f2,fit2));
%                 end
else
fit2 = 0;
end
gfpfit1{t}(dai2,dai) = fit1;
gfpfit2{t}(dai2,dai) = fit2;

%             if(f==conf1(f2) && t==1)
% %             fprintf('t: %i, (%i,%i)-> fit1: %f, fit2: %f\n',t,f,f2,fit1,fit2);
% %             imshowpair(cglow2b,sglow2a);
%                 pause(0.1);
%             end
tots{t}(dai2,dai) = [overs{t}(dai2,dai) areas{t}(dai2,dai)...
peris{t}(dai2,dai) lens{t}(dai2,dai) wids{t}(dai2,dai)...
eyes{t}(dai2,dai) gfpfit1{t}(dai2,dai)...
gfpfit2{t}(dai2,dai)]*weights;
end
end
end

% priority matrix for showing closest matches
priority = cell(1,length(tots));
for t = tmax_lkup:length(tots)
dist = tots{t};
[nr,nc] = size(dist);
pris = zeros(nr,nc);
for i = 1:nr
[~,pris(i,:)] = sort(dist(i,:));
end
priority{t} = pris;
end

% minimize difference measure to identify fishes
% config = cell(1,length(tots));
% for t = 1:length(tots)
%     config{t} = comp_iter2(tots{t});
% end

% using priority matrix, find lowest energy configurations without a
% conflict
for t = tmax_lkup:length(tots)
[nr,~] = size(tots{t});
config{t} = findConfig_v1_3(priority{t},ones(1,nr),tots{t});
end
end

% verify configuration visually
if show == 1
for t = 1:length(BYTIME)-1
conf = config{t};
for dai = 1:length(BYTIME(end).dataAt{t+1})
ff = BYTIME(end).dataAt{t+1}(dai);
conf_ff = BYTIME(end).dataAt{t}(conf(dai));
figure(1);
%             suptitle(['Fish (t' num2str(t+1) ', ' num2str(ff) ') to (t' num2str(t) ', ' num2str(conf_ff) ')']);
subplot(2,2,1);
imshow(imcrop(BYTIME(t).fishes(conf_ff).rounds(3).gfp1,tum_rect),[]);
title({['Fish (t' num2str(t+1) ', ' num2str(ff) ') to (t' num2str(t) ', ' num2str(conf_ff) ')'],['Fish (t' num2str(t) ', ' num2str(conf_ff) ') GFP1']});
subplot(2,2,2);
imshow(imcrop(BYTIME(t).fishes(conf_ff).rounds(3).gfp2,tum_rect),[]);
title(['Fish (t' num2str(t) ', ' num2str(conf_ff) ') GFP2']);
subplot(2,2,3);
imshow(imcrop(BYTIME(t+1).fishes(ff).rounds(3).gfp1,tum_rect),[]);
title(['Fish (t' num2str(t+1) ', ' num2str(ff) ') GFP1']);
subplot(2,2,4);
imshow(imcrop(BYTIME(t+1).fishes(ff).rounds(3).gfp2,tum_rect),[]);
title(['Fish (t' num2str(t+1) ', ' num2str(ff) ') GFP2']);
pause(1);
%             w = waitforbuttonpress;
end
end
end

% the following is backward looking: identifying fish that survived until
% the end and ignoring the rest so the id of the fish is based on the index
% number at the last timepoint
% save configuration arrary at last time point in fish
for ff=1:length(BYTIME(end).fishes)
% see if fish data available
daie = find(BYTIME(end).dataAt{end}==ff);
if(~isempty(daie))
BYTIME(end).fishes(ff).config = zeros(1,length(BYTIME));
for t = 1:length(BYTIME)
dai = daie;
cff = ff;
for k = 1:(length(BYTIME)-t)
kt = length(BYTIME)-k;
dai = config{kt}(dai);
cff = BYTIME(end).dataAt{kt}(dai);
end
BYTIME(end).fishes(ff).config(t) = cff;
end
end
end

% this variable is used later to identify fishes to analyze
BYTIME(end).dataAllTimes = 1:length(BYTIME(end).fishes);

% the following is foward looking: the id of the fish is based on the index
% at the first timepoint
dataUpTo = cell(1,length(BYTIME));

% initialize dataUpTo
for t = 1:length(BYTIME)
dataUpTo{t} = zeros(1,length(BYTIME(end).dataAt{t}));
end
% iterate through time backwards and use config to identify fish
dataCounter = dataUpTo; % binary counter to mark when data has been assigned
for t = 1:length(BYTIME)
tt = length(BYTIME)-t+1;
unused_dai = find(dataCounter{tt}==0);
for dai = unused_dai
dataCounter{tt}(dai) = 1;
% trace down to first timepoint
% assume config points to which indices of dataAt match while
% dataAt contains the actual indices of the fish in BYTIME.fishes
daif = dai;
f = BYTIME(end).dataAt{tt}(daif);
dai_list = zeros(1,tt);
lookup = zeros(1,tt);
dai_list(tt) = daif; % keep track of fish identity over time
lookup(tt) = f;
for k = 1:(tt-1)
kt = tt-k;
daif = config{kt}(daif);
f = BYTIME(end).dataAt{kt}(daif);
dataCounter{kt}(daif) = 1;
dai_list(kt) = daif;
lookup(kt) = f;
end
% mark dataUpTo with identified fish using time 1 label
for ttt = 1:tt
dataUpTo{ttt}(dai_list(ttt)) = f;
end
% save lookup array for fish
BYTIME(1).fishes(f).lookup = lookup;
end
end

% check data counter to make sure everyting is okay
for t = 1:length(BYTIME)
if(~all(dataCounter{t}))
% not all data is accounted off
error('Unaccounted for data at time %i\n',t);
end
end

% save dataUpTO
BYTIME(1).dataUpTo = dataUpTo;
BYTIME(1).dataToUse = dataUpTo{1};

toc;
disp(' ');

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_gender_v1_2(gender_file,specs,show)
%% CreateFishStruct_gender_v1_1
%  Version 1.1
%  Author: Adeyinka Lesi
%  Date: 8/24/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config
%  11/3/19: found an error in the gender assignment code where the gens
%  array below was set to be the length of the total number of time points
%  for all fish rather than the number of timepoints for the individual
%  fish. Since a zero signified a male, this meant that if the number of
%  time points for the individual fish was much smaller than the max, it
%  was possible for a female fish to be assigned male gender since the
%  gender assignment at each time point is averaged to get the overall
%  assignment. This is now fixed. The real problem is female fish that did
%  not survive for a long time (or were not observed for long for another
%  reason) were preferentially set as female

disp('Step 10.6: Assign gender of fish (fixed)');
tic;

% male = 0, female = 1, dead fish = 2
gender_id = dlmread(gender_file);
if(~isempty(specs.selector))
orig_gid = gender_id;
gender_id = zeros(length(specs.selector{1}),length(specs.selector));
for i = 1:length(specs.selector)
gender_id(1:length(specs.selector{i}),i) = orig_gid(specs.selector{i},i);
gender_id(length(specs.selector{i})+1:end,i) = 2;
end
end


for ff = BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(ff).lookup;
gens = zeros(1,length(lookup));
for t = 1:length(lookup)
f = lookup(t);
gens(t) = gender_id(f,t);
BYTIME(t).fishes(f).gender_initial = gens(t);
end
if(any(gens>1))
error('A dead fish got selected for t = %i, ff = %i',...
find(gens>1,1),ff);
end
BYTIME(1).fishes(ff).gender_list = gens;
gender = round(mean(gens));

for t = 1:length(lookup)
f = BYTIME(1).fishes(ff).lookup(t);
BYTIME(t).fishes(f).gender = gender;
end
end

toc;
disp(' ');

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step11_v2_2(specs,show)
%% CreateFishStruct_step11
%  Version 2.2
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 10/20/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% First scale and rotatate fish so that fish at time point 1, 2 and 3 are the same length
% save scaling transform
% use on body, spine and eye and determine landmarks for nonliniar transform
% save non liniar transform
%% Version History
%  1.2: (not keeping changes from 1.1) - switch to using dataUpTo and lookup instead of dataAllTimes and config
%  also changed ctp2trans and imtranfrom to fitgeotrans and imwarp and using
%  bodyEnd instead of spineEnd
%  2.0: changing how transforms are done (trying to avoid non-linear
%  transform)
%  2.1: use eye position as fixed point instead of nose
%  2.2: aligning rgb images as well; need to ensure consistent end points

disp('Step 11: Finding transform functions for making all time images overlay ...' );

tic;

% get width and heigth of images
ww = specs.maximum_width;
hh = specs.maximum_length;
objs(2*length(BYTIME)).body = BYTIME(1).fishes(1).rounds(3).body;
tails(2*length(BYTIME)).edge = [];
for ff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(ff).lookup;
nt = length(lookup);
% retrieve images and get measurments
body_points = zeros(2*nt,4);
edge_pixel_count = zeros(1,2*nt);
for t = 1:nt
% figure out assigning correct fish
cff = lookup(t);

% retrive images from data struc
objs(2*t-1).body  = BYTIME(t).fishes(cff).rounds(3).body;
objs(2*t-1).eye   = BYTIME(t).fishes(cff).rounds(3).eye;
objs(2*t-1).spine = BYTIME(t).fishes(cff).rounds(3).spine;
objs(2*t-1).bf    = BYTIME(t).fishes(cff).rounds(3).bf(:,:,1);

objs(2*t).body  = BYTIME(t).fishes(cff).rounds(3).cbody;
objs(2*t).eye   = BYTIME(t).fishes(cff).rounds(3).ceye;
objs(2*t).spine = BYTIME(t).fishes(cff).rounds(3).cspine;
objs(2*t).bf    = BYTIME(t).fishes(cff).rounds(3).rgb(:,:,1);

% find reference coordinates for all times
objs(2*t-1).end_coor = BYTIME(t).fishes(cff).rounds(3).bodyEnd_coor;
objs(2*t-1).eye_coor = BYTIME(t).fishes(cff).rounds(3).eye_coor;

objs(2*t).end_coor = BYTIME(t).fishes(cff).rounds(3).cbodyEnd_coor;
objs(2*t).eye_coor = BYTIME(t).fishes(cff).rounds(3).ceye_coor;

% input points for length scaling transform (make fish same length on t1, t2 and t3 image)
% nose coor (which is already fixed) and end of spine point
objs(2*t-1).input_points = [objs(2*t-1).eye_coor ; objs(2*t-1).end_coor];
body_points(2*t-1,:) = [objs(2*t-1).eye_coor objs(2*t-1).end_coor];

objs(2*t).input_points = [objs(2*t).eye_coor ; objs(2*t).end_coor];
body_points(2*t,:) = [objs(2*t).eye_coor objs(2*t).end_coor];

% get tail edge image and store
[tails(2*t-1).edge,tails(2*t-1).box] = getTailEdgeImage_v1_0(objs(2*t-1).body,objs(2*t-1).bf,objs(2*t-1).end_coor);
[tails(2*t).edge,tails(2*t).box] = getTailEdgeImage_v1_0(objs(2*t).body,objs(2*t).bf,objs(2*t).end_coor);
edge_pixel_count(2*t-1) = sum(tails(2*t-1).edge(:));
edge_pixel_count(2*t) = sum(tails(2*t).edge(:));
end

% need to align tails - align everything to image with most pixels
[~,imax] = max(edge_pixel_count);
% the bf and rgb images are paired; want the index of the bf image
if(mod(imax,2)==0)
imax = imax - 1;
end
ref1 = tails(imax).edge;
ref2 = tails(imax+1).edge;
[sref2,dr2,dc2] = alignBW_v2_0(ref2,ref1,1,1,10,10);
tails(imax).shifted = ref1;
tails(imax+1).shifted = sref2;
shifts = zeros(2*nt,2);
shifts(imax+1,:) = [dc2 dr2];
% adjust body reference points
objs(imax+1).input_points(2,:) = objs(imax+1).input_points(2,:)-shifts(imax+1,:);
body_points(imax+1,3:4) = body_points(imax+1,3:4)-shifts(imax+1,:);

for t2 = [1:imax-1 imax+2:2*nt]
if(mod(t2,2)==0)
ref = tails(t2-1).shifted;
else
ref = ref1;
end
[tails(t2).shifted,shifts(t2,2),shifts(t2,1)] = alignBW_v2_0(...
tails(t2).edge,ref,1,1,10,10);

% adjust body reference points
objs(t2).input_points(2,:) = objs(t2).input_points(2,:)-shifts(t2,:);
body_points(t2,3:4) = body_points(t2,3:4)-shifts(t2,:);
end

%     % compare images of tails
%     for t = 1:nt
%         figure(2);
%         subplot(1,2,1);
%         imshowpair(tails(2*t-1).edge,tails(2*t).edge);
%         subplot(1,2,2);
%         imshowpair(tails(2*t-1).shifted,tails(2*t).shifted);
%         figure(3);
%         box1 = tails(2*t-1).box; %[objs(2*t-1).input_points(2,:)-[50 60] 70 120];
%         box2 = tails(2*t).box; %[objs(2*t).input_points(2,:)-[50 60] 70 120];
%         tail1 = imcrop(objs(2*t-1).bf,box1);
%         tail2 = imcrop(objs(2*t).bf,box2);
%         subplot(1,2,1);
%         imshow(tail1);
%         hold on;
%         x1 = objs(2*t-1).input_points(2,1)-box1(1)+1;
%         y1 = objs(2*t-1).input_points(2,2)-box1(2)+1;
%         plot(x1,y1,'o');
%         plot(x1+shifts(2*t-1,1),y1+shifts(2*t-1,2),'*');
%         hold off;
%         subplot(1,2,2);
%         imshow(tail2);
%         hold on;
%         x2 = objs(2*t).input_points(2,1)-box2(1)+1;
%         y2 = objs(2*t).input_points(2,2)-box2(2)+1;
%         plot(x2,y2,'o');
%         plot(x2+shifts(2*t,1),y2+shifts(2*t,2),'*');
%         hold off;
%         pause(0.1);
%     end

mean_eye = [mean(body_points(:,1)) mean(body_points(:,2))];
mean_end = [mean(body_points(:,3)) mean(body_points(:,4))];
body_lengths = body_points(:,3)-body_points(:,1);
mean_length = mean_end(1)-mean_eye(1);
std_length = std(body_lengths);
outliers = find(abs(body_lengths-mean_length(1))>2*std_length);

% base points for length scaling transform is eye coor and average
% fish length.
base_points     = [mean_eye ; mean_end];

% It seems from examining the images that fish suffering from melanoma
% might naturally become thinner and less robust looking over time - it
% might be inproper to try and overlay body images over time in this
% case because the assumption that the body size is approximately the
% same over time is violated. It might be best to rely mostly on the
% spine, eye, nose and tail positions (The spine identification as been
% improved to 1) avoid getting messed up due to pigmented melanoma
% cells and 2) to extend further into the body of the fish - the belly
% has been removed from the bw spine image 3) the spine is also more or
% less represented as a smooth curve

% return input and base points from fish body boundaries
%     [x_coor,y_updown,y_coor_base] = return_front_coor_T_mod2(objs(1:nt));
% return input and base points along spine in tail
[sx_coor,sy_coor,sx_coor_base,sy_coor_base] = return_spine_coor_T_mod2(objs(1:2*nt),8);
% other points to help with transformation - these points define
% alignment in the y direction - all the points are spaced the same
% amount so that mostly translation and not stretching occurs in
% the y direction.
up_spacing = 150;
down_spacing = 200;
center_points_x = zeros(2*nt,length(base_points(:,1))+length(sx_coor_base));
center_points_y = center_points_x;
for t2 = 1:2*nt
center_points_x(t2,:) = [objs(t2).input_points(:,1)' sx_coor(t2,:)];
center_points_y(t2,:) = [objs(t2).input_points(:,2)' sy_coor(t2,:)];
end
center_base_x = [base_points(:,1)' sx_coor_base];
center_base_y = [base_points(:,2)' sy_coor_base];
y_up = center_points_y-up_spacing;
y_down = center_points_y+down_spacing;
y_up_base = center_base_y-up_spacing;
y_down_base = center_base_y+down_spacing;

for t2 = 1:2*nt
objs(t2).input_points_T = [[center_points_x(t2,:)',y_up(t2,:)'];...
[center_points_x(t2,:)',center_points_y(t2,:)'];[center_points_x(t2,:)',y_down(t2,:)']];
end

base_points_T = [[center_base_x',y_up_base'];...
[center_base_x',center_base_y'];[center_base_x',y_down_base']];

for t = 1:nt
% figure out assigning correct fish
cff = lookup(t);

% calculate image tranformation (nonliniar)
TFORM_T = fitgeotrans(objs(2*t-1).input_points_T,base_points_T,'polynomial',2);
TFORM_T2 = fitgeotrans(objs(2*t).input_points_T,base_points_T,'polynomial',2);

% save transform functions in data struct.
BYTIME(t).fishes(cff).transforms.timeTF_bf = TFORM_T;
BYTIME(t).fishes(cff).transforms.timeTF_rgb = TFORM_T2;

% imref2d object for imwarp
RA = imref2d(size(objs(2*t-1).body),[0.5 hh+0.5],[0.5 ww+0.5]);
% transform bw images (scaling and rotation)
objs(2*t-1).body_T   = logical(imwarp(objs(2*t-1).body,TFORM_T,'OutputView',RA,'FillValues',0));
objs(2*t-1).spine_T  = logical(imwarp(objs(2*t-1).spine,TFORM_T,'OutputView',RA,'FillValues',0));
objs(2*t-1).eye_T    = logical(imwarp(objs(2*t-1).eye,TFORM_T,'OutputView',RA,'FillValues',0));

objs(2*t).body_T   = logical(imwarp(objs(2*t).body,TFORM_T2,'OutputView',RA,'FillValues',0));
objs(2*t).spine_T  = logical(imwarp(objs(2*t).spine,TFORM_T2,'OutputView',RA,'FillValues',0));
objs(2*t).eye_T    = logical(imwarp(objs(2*t).eye,TFORM_T2,'OutputView',RA,'FillValues',0));

% calculate distortion
BYTIME(t).fishes(cff).transforms.timeTF_bf_distortion = sum(objs(2*t-1).body_T(:))/sum(objs(2*t-1).body(:));
BYTIME(t).fishes(cff).transforms.timeTF_rgb_distortion = sum(objs(2*t).body_T(:))/sum(objs(2*t).body(:));
end

if show==1       
sum_body = zeros(size(objs(1).spine));
sum_body2 = zeros(size(objs(1).spine_T));
for t2 = 1:2*nt
% outlines for plot 
objs(t2).BB1 =  bwboundaries(objs(t2).body);
sum_body = sum_body + objs(t2).spine-objs(t2).eye;
objs(t2).B1 =  bwboundaries(objs(t2).body_T);
sum_body2 = sum_body2 + objs(t2).spine_T-objs(t2).eye_T;
end

figure(1)
subplot(2,1,1)
imshow(sum_body,[])
hold on
title(['All ' num2str(nt) ' time points, before transform. id: ' num2str(ff) ...
', mean length: ' num2str(round(mean_length)) ' ± ' ...
num2str(round(std_length))],'FontSize',10)
cols = 'cbkyg';
for t2 = 1:2*nt
for k=1:length(objs(t2).BB1)
b=objs(t2).BB1{k};
plot(b(:,2),b(:,1),['-' cols(mod(t2-1,5)+1)])
end
plot(objs(t2).input_points_T(:,1),objs(t2).input_points_T(:,2),'+r','MarkerSize',15);
end
plot(base_points(:,1),base_points(:,2),'og','MarkerSize',15)
%         legend('t = 1, outline','t = 2, outline','t = 3, outline','Input points','Base point','Location','NorthEastOutside')
hold off;

subplot(2,1,2)
imshow(sum_body2,[])
hold on
title(['All ' num2str(nt) ' time points, after transform. id: ' num2str(ff) ],'FontSize',10)
for t2 = 1:2*nt
for k=1:length(objs(t2).B1)
b=objs(t2).B1{k};
plot(b(:,2),b(:,1),['-' cols(mod(t2-1,5)+1)])
end
end
plot(base_points(:,1),base_points(:,2),'og','MarkerSize',15)
plot(base_points_T(:,1),base_points_T(:,2),'og','MarkerSize',15)
hold off;
if(~isempty(outliers))
disp(['Outliers at id=' num2str(ff) ', t=' num2str(outliers')]);
end
pause(0.1);
end

end
toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step12_v2_1(specs,show)
%% CreateFishStruct_step12
%  Version 1.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 9/12/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History:
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config
%  switched from imtransform to imwarp
%  2.0: changing transforms used - only a 'slight' polynomial transform
%  2.1: processing rgb images

%% Use transforms found in previous step. Save transformed images in rounds(4)
disp('Step 12: Using transform functions found in previous step for making all three time point images overlay ... Saving transformed images in rounds(4)' );

tic;

% get width and heigth of images
ww = specs.maximum_width;
hh = specs.maximum_length;

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t= 1:length(lookup)
% figure out assigning correct fish
ff = lookup(t);

% retrive images from data struc
body  = BYTIME(t).fishes(ff).rounds(3).body;
eye  = BYTIME(t).fishes(ff).rounds(3).eye;
spine  = BYTIME(t).fishes(ff).rounds(3).spine;
bf  = BYTIME(t).fishes(ff).rounds(3).bf;
gfp1  = BYTIME(t).fishes(ff).rounds(3).gfp1;
gfp2  = BYTIME(t).fishes(ff).rounds(3).gfp2;
rfp1  = BYTIME(t).fishes(ff).rounds(3).rfp1;
rfp2  = BYTIME(t).fishes(ff).rounds(3).rfp2;

cbody = BYTIME(t).fishes(ff).rounds(3).cbody;
ceye = BYTIME(t).fishes(ff).rounds(3).ceye;
cspine = BYTIME(t).fishes(ff).rounds(3).cspine;
rgb = BYTIME(t).fishes(ff).rounds(3).rgb;

% retriave transform functions from data struct.
timeTF_bf = BYTIME(t).fishes(ff).transforms.timeTF_bf;
timeTF_rgb = BYTIME(t).fishes(ff).transforms.timeTF_rgb;

% get imref2d object
RA = imref2d(size(body),[0.5 hh+0.5],[0.5 ww+0.5]);
bf_fill = round(mean(bf(bf>mean(bf(:)))));

body_TT   = logical(imwarp(body,timeTF_bf,'OutputView',RA,'FillValues',0));
spine_TT  = logical(imwarp(spine,timeTF_bf,'OutputView',RA,'FillValues',0));
eye_TT    = logical(imwarp(eye,timeTF_bf,'OutputView',RA,'FillValues',0));
bf_TT     = imwarp(bf,timeTF_bf,'OutputView',RA,'FillValues',bf_fill);
gfp1_TT   = imwarp(gfp1,timeTF_bf,'OutputView',RA,'FillValues',0);
gfp2_TT   = imwarp(gfp2,timeTF_bf,'OutputView',RA,'FillValues',0);
rfp1_TT   = imwarp(rfp1,timeTF_bf,'OutputView',RA,'FillValues',0);
rfp2_TT   = imwarp(rfp2,timeTF_bf,'OutputView',RA,'FillValues',0);

cbody_TT   = logical(imwarp(cbody,timeTF_rgb,'OutputView',RA,'FillValues',0));
cspine_TT  = logical(imwarp(cspine,timeTF_rgb,'OutputView',RA,'FillValues',0));
ceye_TT    = logical(imwarp(ceye,timeTF_rgb,'OutputView',RA,'FillValues',0));
rgb_TT = rgb;
for dim = 1:size(rgb,3)
temp = rgb(:,:,dim);
temp_fill = round(mean(temp(temp>mean(temp(:)))));
tempTT = imwarp(temp,timeTF_rgb,'OutputView',RA,'FillValues',temp_fill);
rgb_TT(:,:,dim) = tempTT;
end

% save images in data struc in rounds(4)
BYTIME(t).fishes(ff).rounds(4).body  = body_TT;
BYTIME(t).fishes(ff).rounds(4).eye   = eye_TT;
BYTIME(t).fishes(ff).rounds(4).spine = spine_TT;
BYTIME(t).fishes(ff).rounds(4).bf    = bf_TT;
BYTIME(t).fishes(ff).rounds(4).gfp1  = gfp1_TT;
BYTIME(t).fishes(ff).rounds(4).gfp2  = gfp2_TT;
BYTIME(t).fishes(ff).rounds(4).rfp1  = rfp1_TT;
BYTIME(t).fishes(ff).rounds(4).rfp2  = rfp2_TT;

BYTIME(t).fishes(ff).rounds(4).cbody  = cbody_TT;
BYTIME(t).fishes(ff).rounds(4).ceye   = ceye_TT;
BYTIME(t).fishes(ff).rounds(4).cspine = cspine_TT;
BYTIME(t).fishes(ff).rounds(4).rgb    = rgb_TT;

% visually inspect transformed images
if show==1
figure(1)
subplot(2,2,1)
imshow(body_TT+spine_TT+eye_TT,[])
title(['(id=' num2str(fff) ',f=' num2str(ff) ',t=' num2str(t) ')'])
subplot(2,2,2)
imshow(bf_TT,[])
title('non-linear transforms')
subplot(2,2,3)
imshow(rfp1_TT+rfp2_TT,[])
title('RFP')
subplot(2,2,4)
imshow(gfp1_TT+gfp2_TT,[])  
title('GFP')  
figure(2);
subplot(2,1,1)
imshow(cbody_TT+cspine_TT+ceye_TT,[]);
title(['(id=' num2str(fff) ',f=' num2str(ff) ',t=' num2str(t) ')'])
subplot(2,1,2)
imshow(rgb_TT,[]);
title('non-linear transforms')
pause(0.1);
end

end
end

toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step13_v2_9(specs,name,show)
%% CreateFishStruct_step13_v2_8
%  Version 2.8
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 9/21/19
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Subtract autoflourescence using RFP image and segment gfp expressing regions in fish 
%%  Version History:
%  1.1: using varying threshold detection (return_side_mod2)
%  1.1: not using varying threshold (return_side_mod); saving image to file
%  1.3: using return_side
%  1.4: return_side_mod3
%  1.5: remove tail and head artifacts
%  1.6: switch to using dataUpTo and lookup instead of dataAllTimes and config
%  1.7: change filter for tail artifacts, to remove less; and how much of
%  body border is removed; also, changed return_side_mod; added rfp mask to
%  make sure bright AF is removed
%  2.0: adjusting tumor recognition method
%  2.1: dilating rbw images, saving bw_pig
%  2.2: changes to limit false positives, and identify all tumor area, may
%  use rgb image info
%  2.3: return_side_mod9: using low threshold as well
%  2.4: subtractiong AF from region with strong AF signal separately from
%  rest (with a higher multiplier). Also Raising removal threshold
%  2.5: return_side_mod10: using reference for fluorescent images
%  2.6: changes for Silja data to improve thresholding and recognition of
%  pigmentation
%  2.7: took some day 1 data - turns out the pigmentation detection
%  sometimes picks up the backbone if the bf image isn't bright or picks up
%  wounds (especially at 0DPT). Will set bf_black to blank at first two
%  time points
%  2.8: using return_side_mod11
%  2.9: using return_side_mod13 and adding brightness criteria for removing
%  AF through the region method

disp('Step 13: Subtracting autoflourescence, segmenting images and saving in rounds(4)...' );

tic;

for fff= BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t= 1:length(lookup)
% figure out assigning correct fish
ff = lookup(t);

% retrive images from data struct
body  = BYTIME(t).fishes(ff).rounds(4).body;
eye  = BYTIME(t).fishes(ff).rounds(4).eye;
gfp1  = BYTIME(t).fishes(ff).rounds(4).gfp1;
gfp2  = BYTIME(t).fishes(ff).rounds(4).gfp2;
rfp1  = BYTIME(t).fishes(ff).rounds(4).rfp1;
rfp2  = BYTIME(t).fishes(ff).rounds(4).rfp2;
bf = BYTIME(t).fishes(ff).rounds(4).bf(:,:,1);
thinbody = bwmorph(body,'thin',20);
thinbody = bwmorph(thinbody,'spur',20);

% find boundary points around eye for plotting
E = bwboundaries(eye);

% erode eye as to not remove metastasis around eye
eye = bwmorph(eye,'erode',5);

% reference to normalize images
ref_int = specs.reference_threshold*65535;

% make sure bright autofluoresence pixels are not included (usually
% the belly of the fish)
% 1) use return_side on rfp images
[rbw1a,~] = return_side_mod13(rfp1,body,eye,specs,100,2,max(ref_int,max(rfp1(:))));
[rbw2a,~] = return_side_mod13(rfp2,body,eye,specs,100,2,max(ref_int,max(rfp2(:))));
% 2) use thresholding on rfp images
% left image
%         rbw1b = binarizeByPercentile_v1_0(rfp1,body,0.999);
%         % clean up
%         rbw1b = bwmorph(rbw1b,'dilate',2);
%         rbw1b = bwmorph(rbw1b,'erode',3);
%         rbw1b = bwmorph(rbw1b,'dilate',2);
%         rbw1b = imfill(rbw1b,'holes');
%         % right image
%         rbw2b = binarizeByPercentile_v1_0(rfp2,body,0.999);
%         % clean up
%         rbw2b = bwmorph(rbw2b,'dilate',2);
%         rbw2b = bwmorph(rbw2b,'erode',3);
%         rbw2b = bwmorph(rbw2b,'dilate',2);
%         rbw2b = imfill(rbw2b,'holes');

% combine two rbw images
rbw1 = rbw1a;
%         rbw1(rbw1b==1) = 1;
rbw2 = rbw2a;
%         rbw2(rbw2b==1) = 1;
% dilate rbw images
rbw1 = bwmorph(rbw1,'dilate',2);
rbw2 = bwmorph(rbw2,'dilate',2);
rbw = rbw1;
rbw(rbw2==1) = 1;

%subtract autofluorescence signal        
af_factor = max(0,specs.rfp_gfp_ratio-1);
af1 = rfp1;
af1(rbw1==0) = 0;
af2 = rfp2;
af2(rbw2==0) = 0;

im1 = gfp1-rfp1-af1*af_factor;
im2 = gfp2-rfp2-af2*af_factor;

% want to compile real tumor regions from previous times
% due to noise, will remove anything close to edges
past_tumors = false(size(im1));
for tt = 1:t-1
past_tumors(BYTIME(tt).fishes(lookup(tt)).rounds(4).bw) = 1;
end
past_tumors(thinbody==0) = 0;

% function 'return_side()' tries to identify actual tumor regions
% from the signal
[bw1,candidates1] = return_side_mod13(im1,body,eye,specs,1000,4,max(im1(:)),past_tumors);
[bw2,candidates2] = return_side_mod13(im2,body,eye,specs,1000,4,max(im2(:)),past_tumors);


% check to see if there are extra bright regions that shouldn't be
% removed by removeOverlappingRegions_v1_0
im1_high_thresh = imbinarize(im1,specs.bw_threshold_high);
im2_high_thresh = imbinarize(im2,specs.bw_threshold_high);
im_high_thresh = im1_high_thresh|im2_high_thresh;

% keep track of pixels removed
removed_count = zeros(1,5);
% remove autofluorescence and merge left and right segmented images
bw1c = removeOverlappingRegions_v1_0(bw1,rbw-im_high_thresh,0.5);
rem_bw1 = bw1-bw1c;
bw2c = removeOverlappingRegions_v1_0(bw2,rbw-im_high_thresh,0.5);
rem_bw2 = bw2-bw2c;
% remove regions found to be autofluoresence in other image
bw1c = removeOverlappingRegions_v1_0(bw1c,rem_bw2-im_high_thresh,0.75);
bw2c = removeOverlappingRegions_v1_0(bw2c,rem_bw1-im_high_thresh,0.75);
rem_bw1 = bw1-bw1c;
removed_count(1) = sum(rem_bw1(:));
rem_bw2 = bw2-bw2c;
removed_count(2) = sum(rem_bw2(:));

% tumors with high pigmentation are hard to detect entirely since
% signal is obscured - will look for black pixels in the bf image
% to find the eye and any pigmented cells. The eye serves as a good
% measure for how black the tumor will be so we will derive a
% threshold from the eye color
if(t>=specs.pigmentation_use_index) % only use at later time points
bf_eye = double(bf(eye==1))/65535;
mean_bf_eye = mean(bf_eye);
std_bf_eye = std(bf_eye);
black_thresh = max(mean_bf_eye+4*std_bf_eye,specs.min_pigmentation_threshold);
bf_black = 1-imbinarize(bf,black_thresh);
% remove eye area and region behind eye since those are dark and
% will show up
bf_black(bwmorph(eye,'dilate',10)==1) = 0; % remove eye
bf_black(:,1:specs.eye_coordinate(1)+specs.eye_diameter*6) = 0; % having some issue with dark areas near eye - so removing a lot
% clean up binarized image (want to be aggressive with this to
% remove false positives)
bf_black = bwmorph(bf_black,'dilate',5);
bf_black = bwmorph(bf_black,'erode',5);
bf_black = bwareaopen(bf_black,20);
bf_black = imfill(bf_black,'holes');
else
bf_black = false(size(bf));
end

% combine bw images
bw = bw1c;
bw(bw2c==1)=1;
% as an additional check, will only use regions that intersect with
% detected signal - the pigmentation will only serve to extend
% known tumors rather than create new ones;
bw_pig = bwlabel(bf_black);
bw_pig_valid = false(size(bw_pig));
nregs = max(bw_pig(:));
for n = 1:nregs
if(any(bw(bw_pig==n)))
bw(bw_pig==n)=1;
bw_pig_valid(bw_pig==n)=1;
end
end

% remove small objects
bw = bwareaopen(bw,5);
bw = imfill(bw,'holes');

% remove tail section and mouth from bw image
bw_orig = bw;
[~, tail_col] = find(thinbody,10,'last');
mean_tail_col = round(mean(tail_col));
tail_cut_col = mean_tail_col-10;
eye_coor = round(0.5*(BYTIME(t).fishes(ff).rounds(2).eye1_coor...
+ BYTIME(t).fishes(ff).rounds(2).eye2_coor));
% remove clusters near tail end
center_tail = false(size(bw));
center_tail(:,tail_cut_col:end) = 1;
bw = removeOverlappingRegions_v1_0(bw,center_tail,0.33);
rem_bw3 = bw_orig-bw;
removed_count(3) = sum(rem_bw3(:));
bw_orig = bw;
% remove clusters near eye and mouth
mouth = false(size(bw));
mouth(eye_coor(2):end,1:eye_coor(1)) = 1;
for r = 1:eye_coor(2)-1
c = r+eye_coor(1)-eye_coor(2);
if(c>0)
mouth(r,1:c) = 1;
end
end
bw = removeOverlappingRegions_v1_0(bw,mouth,0.33);
rem_bw4 = bw_orig-bw;
removed_count(4) = sum(rem_bw4(:));
bw_orig = bw;
% screen outer tail regionn
% a cutout of the region around the tail end of the fish
outer_tail = body;
outer_tail(thinbody==1) = 0;
outer_tail_cut_col = mean_tail_col - 20;
outer_tail(:,1:outer_tail_cut_col) = 0;
bw = removeOverlappingRegions_v1_0(bw,outer_tail,0.25);
rem_bw5 = bw_orig - bw;
removed_count(5) = sum(rem_bw5(:));

bw_dif = rem_bw1+rem_bw2+rem_bw3+rem_bw4+rem_bw5;

% take brigthtest pixel from the left and right pure signal images
% and make merged image
im = Merge(im1,im2);

% save merged pure signal and merged segmented image in rounds(4) 
BYTIME(t).fishes(ff).rounds(4).bw = logical(bw);
BYTIME(t).fishes(ff).rounds(4).pure = uint16(im);
BYTIME(t).fishes(ff).rounds(4).bw_candidates1 = candidates1;
BYTIME(t).fishes(ff).rounds(4).bw_candidates2 = candidates2;
BYTIME(t).fishes(ff).rounds(4).omitted_pixel_counts = removed_count;
BYTIME(t).fishes(ff).rounds(4).bw_pigmented = logical(bw_pig_valid);

% Show images with autoflourescence subtracted                      
B = bwboundaries(body);
C = bwboundaries(bw);
D = bwboundaries(outer_tail);
F = bwboundaries(bw_dif);
G1 = bwboundaries(rbw1);
G2 = bwboundaries(rbw2);
M = bwboundaries(mouth);

fig1 = figure(1);
if show == 0
fig1.Visible = 'off';
end
subplot(1,1,1);
imshow(double(im)/double(max(im(:))));
hold on
title(...
['GFP signal, merged, minus autoflourescence. id: ' num2str(fff) ', f: ' ...
num2str(ff) ',  t: ' num2str(t) ', Area Omitted: ' ...
num2str(removed_count)],'FontSize',13)
plot([tail_cut_col tail_cut_col],[1 size(im,1)],'y-','Linewidth',1.25);
%         plot([1 eye_coor(1)],[eye_coor(2) eye_coor(2)],'y-','Linewidth',1.25);
%         plot([eye_coor(1) eye_coor(1)],[eye_coor(2) size(im,1)],'y-','Linewidth',1.25);
b=B{1};
plot(b(:,2),b(:,1),'b-','LineWidth',5)
b=E{1};
plot(b(:,2),b(:,1),'r-','LineWidth',3)
for c=1:length(C)
b=C{c};
plot(b(:,2),b(:,1),'g-')
end
b = D{1};
plot(b(:,2),b(:,1),'y-','LineWidth',1.5)
for c=1:length(F)
b=F{c};
plot(b(:,2),b(:,1),'r--')
end
for c=1:length(G1)
b=G1{c};
plot(b(:,2),b(:,1),'y--')
end
for c=1:length(G2)
b=G2{c};
plot(b(:,2),b(:,1),'y--')
end
for c=1:length(M)
b=M{c};
plot(b(:,2),b(:,1),'y-','LineWidth',1.5)
end

if(~exist(name,'dir'))
mkdir(name);
end
img_name = ['id' num2str(fff) '_t' num2str(t) '_' BYTIME(t).fishes(ff).fishID];
saveas(fig1,[name '/' img_name],'jpg');

if(show == 1)
pause(0.1);
end
clf;
end
end

toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step14_v1_1(show)
%% CreateFishStruct_step14
%  Version 1.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 8/24/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Mask day 1 using day 7 bw image and mask day 7 bw image using day 14 bw image
%we assume that object not segmented a week later where noise when they
%where first detected. Save in rounds(4) (NB thus overwriting previous image)
%% Version History
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config

disp('Step 14: Removing objects which are no longer there in the next timepoint...' );

tic;

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
% retrieve images and get measurments
objs(length(lookup)).bw_t1 = [];
for t = 1:length(lookup)
% figure out assigning correct fish
ff = lookup(t);

% retrieve images from data struc
objs(t).bw_t1  = BYTIME(t).fishes(ff).rounds(4).bw;
objs(t).im_t1  = BYTIME(t).fishes(ff).rounds(4).pure;
objs(t).mask = bwmorph(objs(t).bw_t1,'dilate',1);
end

for t = 1:length(lookup)-1
% mask_t2 should not have objects which are not present in mask_t3
objs(t).mask(objs(t+1).mask==0)=0;
objs(t).clean_bw_t1 = objs(t).bw_t1;
% clean up the past
objs(t).clean_bw_t1(objs(t+1).mask==0)=0;
end
objs(end).clean_bw_t1 = objs(end).bw_t1;

% show final segmentation result on top of pure=gfp-rfp images
if show==1
npic = length(lookup);
figure(1)
for t = 1:npic
C1 = bwboundaries(objs(t).clean_bw_t1);

subplot(npic,1,t)
imshow(objs(t).im_t1.*1)
hold on

% figure out assigning correct fish
ff = lookup(t);

title(['Green line: segmentation. Fish #: ' num2str(ff)  ',  Time: ' num2str(t)],'FontSize',13)
for k=1:length(C1)
b=C1{k};
plot(b(:,2),b(:,1),'-g')
end
end
pause(1);
end

% save cleaned up segmentation from t=1 and t=2 images in rounds(4)
for t = 1:length(lookup)
% figure out assigning correct fish
ff = lookup(t);
BYTIME(t).fishes(ff).rounds(4).clean_bw = logical(objs(t).clean_bw_t1);
end
end


toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step15_v1_1(specs,show)
%% CreateFishStruct_step15
%  Version 1.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 8/24/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%%  Construct DistMatrix form BW segmented images (pairwise distance
% matrix holding the shortest distances between all seperate regions with
% 1's. (distMatrix is a symmetric matrix))
% Distances are shortest possible distances between points on the
% boundary of both objects
% use this matrix to make labeled segmented image where each cluster of
% close objects have same value 1,2,3 ect and the rest is 0
%% Version History:
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config

disp('Step 15: Finding pairwise distance matrix for all connected regions in segmented images and use this to cluster close objects...' );

tic;

% threshold (in pixels) for clustering close objects
threshold = specs.clustering_threshold;

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t = 1:length(lookup)
% figure out assigning correct fish
ff = lookup(t);

% retrive images from data struc
bw = BYTIME(t).fishes(ff).rounds(4).bw;

% function 'return_dist_matrix' finds matrix holding the shortest distances between all seperate regions with
% 1's. Distances are shortest possible distances between points on the
% boundary of both objects. If there is only one object in fish (only
% primary) distmatrix is empty []
distMatrix = return_dist_matrix(bw);        

L_clustered = return_clustered_L_mod(bw,distMatrix,threshold,t);

% save matrix and L_clustered in data struct
BYTIME(t).fishes(ff).rounds(4).distMatrix = distMatrix;
BYTIME(t).fishes(ff).rounds(4).L = L_clustered;

if show==1
figure(1)
ncols = ceil(length(lookup)/4);
subplot(ceil(length(lookup)/ncols),ncols,t)
imshow(double(L_clustered)/double(max(L_clustered(:)))*63,colormap('jet'));
%             hold on
title(['id: ' num2str(fff)  ', t: ' num2str(t)],'FontSize',10);
if(t==length(lookup))
pause(1);
end
end
end
end

toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step16_v1_5(specs,name,show)
%% Separating merged objects
%% CreateFishStruct_step16
%  Version 1.5
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 10/17/19
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  1.1: using return_L_divided_mod4
%  1.2: switch to using dataUpTo and lookup instead of dataAllTimes and config
%  1.3: save images to file
%  1.4: different Day1 clustering threshold
%  1.5: using return_L_divided_mod7 (better segmentation)

disp('Step 16: Separating objects which has merged due to melanoma growth since last timepoint...' );
tic;


objs(length(BYTIME)) = struct();
for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
nt = length(lookup);
for t = 1:nt
% figure out assigning correct fish
ff = lookup(t);
% retrive images from data struc
objs(t).L_t1 = BYTIME(t).fishes(ff).rounds(4).L;
if(t > 1)
objs(t).threshold = specs.clustering_threshold;
else
objs(t).threshold = specs.clustering_threshold_day1;
end
objs(t).fishID = ['id' num2str(fff) '_' BYTIME(t).fishes(ff).fishID];
end

% function return_L_divided (uses the functions: return_met_events, return_clustered_C, return_object_voronoi_points and return_labeled_voronoi)
% The function seperates objects in labeled matrix L_t which started out as
% seperated object in L_t-1 but has merged over time. When two objects have fully merged (no 0's between then) seperation line is found using
% dilation of ojects in the previous timepoints and using the
% boundary points of these objects to determine voronoi regions
% and assign these regions to either of the two objects
% (see more information in the individual functions listed above)
[L_divideds,V_merges,cross_index,cross_links] = return_L_divided_mod7(objs(1:nt));


fig1 = figure(1);
if show == 0 
fig1.Visible = 'off';
end
fig1.Position = [962 42 958 954];
ncols = ceil(nt/4);

for t = 1:nt
subplot(ceil(nt/ncols),ncols,t)
imshow(L_divideds{t}*(length(hot)-1)/max(L_divideds{t}(:)),hot)
hold on
title(['id: ' num2str(fff)  ', t: ' num2str(t) ', count: '...
num2str(max(L_divideds{t}(:)))],'FontSize',10);
end

hold off;
img_name = ['id' num2str(fff) '_' BYTIME(1).fishes(fff).batch '_' ...
BYTIME(1).fishes(fff).implantSize '_' BYTIME(1).fishes(fff).location];
saveas(fig1,[name '/' img_name],'jpg');
if show==1
pause(0.1);
end
clf;

for t = 1:nt
% figure out assigning correct fish
ff = lookup(t);
% save labeled segmented image in data struct.
BYTIME(t).fishes(ff).rounds(4).L_divided = int8(L_divideds{t});
BYTIME(t).fishes(ff).cross_index = cross_index{t};
end

end


toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_eliminations_v1_1(elimination_file,show)
%% CreateFishStruct_eliminations_v1_1
%  Version 1.1
%  Author: Adeyinka Lesi
%  Date: 8/24/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History:
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config

disp('Step 16.5: Create a list of fish to be analyzed by eliminating selected fish');
tic;

elimCounter = zeros(size(BYTIME(1).dataUpTo{1}));
if(~isempty(elimination_file))
% format: time_index, fish_number
elims = dlmread(elimination_file);
for i = 1:size(elims,1)
for ff = BYTIME(1).dataUpTo{1}
if(BYTIME(1).fishes(ff).lookup(elims(i,1))==elims(i,2))
% this fish is to be eliminated
elimCounter(ff) = 1;
end
end
end
end

BYTIME(1).dataToEliminate = BYTIME(1).dataUpTo{1}(elimCounter==1);
BYTIME(1).dataToUse = BYTIME(1).dataUpTo{1}(elimCounter==0);
toc;
disp(' ');

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ aveFish] = CreateFishStruct_step17_v1_1(specs,show)
%% CreateFishStruct_step17
%  Version 1.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 8/24/16
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Transform to standard fish shape
% First transform all fish to have same length
% the determine average fish shape and transform all images to this
% shape. If the boundary of a fish bulges out over the standard fish shape
% close to the injection site we do not use standard shape as base coor in
% this region - this way we avoid transforming a bulging tumor smaller 
%% Version History:
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config

disp('Step 17: Finding average fish shape of entire data set and transforming to this shape...' );
tic;

% find average fish length
FL=zeros(1,length(BYTIME(1).dataToUse)*length(BYTIME));% vector for storing fish lengths
cc=0; % index for fish length vector

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t=1:length(lookup)
ff = lookup(t);
body  = BYTIME(t).fishes(ff).rounds(4).body;        
cc=cc+1;
FL(cc)= return_fish_length(body);
end
end

% find average fish length
ave_fish_length = round(mean(FL(1:cc)));
% nose_coor_fix already defined
nose_coor_fix = specs.nose_coordinate; % pick coordinates very close to avearge 'nose' location
end_spine_coor_fix = [ave_fish_length, nose_coor_fix(2)-20];

% get width and heigth of images
ww = specs.maximum_width;
hh = specs.maximum_length;

% initialize sum images
bodySum  = zeros(size(BYTIME(end).fishes(1).rounds(4).body));
eyeSum   = zeros(size(BYTIME(end).fishes(1).rounds(4).body));
spineSum = zeros(size(BYTIME(end).fishes(1).rounds(4).body));

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t=1:length(lookup)
ff = lookup(t);
% retrive images from data struc
body   = BYTIME(t).fishes(ff).rounds(4).body;
spine  = BYTIME(t).fishes(ff).rounds(4).spine;
eye    = BYTIME(t).fishes(ff).rounds(4).eye;
bf     = BYTIME(t).fishes(ff).rounds(4).bf;

% find coordinates of end of spine
[y,x]=find(spine,20,'last');
end_spine_coor = [round(mean(x)),round(mean(y))];

% calculate image tranformation (nonliniar)
TFORM = cp2tform([nose_coor_fix;end_spine_coor], [nose_coor_fix;end_spine_coor_fix],'nonreflective similarity');

% save liniar transform function in data struct.
BYTIME(t).fishes(ff).transforms.lengthScaleTF2 = TFORM;

% transform all fish bodies to be the same length
body_T  = logical(imtransform(body,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
eye_T   = logical(imtransform(eye,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
spine_T = logical(imtransform(spine,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
bf_T    = uint16(imtransform(bf,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));

% add all length-rescaled fish bodies, eye and spine in bodySum,
% eyeSum and spineSum images
bodySum  = bodySum + body_T - eye_T;
eyeSum   = eyeSum + eye_T;
spineSum = spineSum + spine_T;

if show==1
figure(1)
imshow(bodySum,[])
hold on
plot(end_spine_coor(1),end_spine_coor(2),'b*','MarkerSize',20)
plot(end_spine_coor_fix(1),end_spine_coor_fix(2),'r*','MarkerSize',20)
legend('Input points for lengthScaleTF2','Base points for lengthScaleTF2')
plot(nose_coor_fix(1),nose_coor_fix(2),'b*','MarkerSize',20)
plot(nose_coor_fix(1),nose_coor_fix(2),'r*','MarkerSize',20)
title('Finding average fish shape and average eye shape and eye location','FontSize',15);
pause(1);
end
end
end

% Using the image bodySum find region where 50% or more of the rescaled fish bodies in the dataset overlap
% this will be the standard fish shape which we map all bodies to
bodySum = bodySum./max(bodySum(:));
eyeSum = eyeSum./max(eyeSum(:));
spineSum = spineSum./max(spineSum(:));
aveFish = smoothBW_mod(imbinarize(bodySum,0.5),30);
aveEye = smoothBW_mod(imbinarize(eyeSum,0.5),10);
aveSpine = smoothBW_mod(imbinarize(spineSum,0.9),10);

% find center of mass of average eye
r = regionprops(aveEye,'Centroid');
aveEyeCoor = [r.Centroid];

if show==1 % plot stadard fish shape
A = bwboundaries(aveFish);
E = bwboundaries(aveEye);
S = bwboundaries(aveSpine);
hold on
b=A{1};
plot(b(:,2),b(:,1),'-y','LineWidth',2)
b=E{1};
plot(b(:,2),b(:,1),'-y','LineWidth',2)
plot(aveEyeCoor(1),aveEyeCoor(2),'oy','LineWidth',2)
hold on
b=S{1};
plot(b(:,2),b(:,1),'-y','LineWidth',2)
pause(0.1);
end

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t=1:length(lookup)
ff = lookup(t);
% retrive images from data struct.
body   = BYTIME(t).fishes(ff).rounds(4).body;
spine   = BYTIME(t).fishes(ff).rounds(4).spine;
eye   = BYTIME(t).fishes(ff).rounds(4).eye;
bf   = BYTIME(t).fishes(ff).rounds(4).bf;

% retrive liniar transform function in data struct.
TFORM1 = BYTIME(t).fishes(ff).transforms.lengthScaleTF2;

% transform all fish bodies to be the same length
body_T = logical(imtransform(body,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
spine_T = logical(imtransform(spine,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
eye_T = logical(imtransform(eye,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
bf_T = imtransform(bf,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);

% find center of mass of eye
r = regionprops(eye_T,'Centroid');
eye_T_Coor = [r.Centroid];

% find input points for transform from fish body and spine in tail
[input_points,base_points] = return_front_coor_Std(body_T,aveFish);
[x_coor,y_coor,y_coor_base] = return_spine_coor_Std(spine_T,aveSpine); % uses function return_Y_coor()

% concatenate points
input_points = [input_points ;[x_coor; y_coor]'     ;eye_T_Coor];
base_points  = [base_points  ;[x_coor; y_coor_base]';aveEyeCoor];

% calculate image tranformation (nonliniar)
TFORM2 = cp2tform(input_points, base_points,'polynomial',2);

% save non liniar transform function in data struct.
BYTIME(t).fishes(ff).transforms.standardFishTF = TFORM2;

if show==1
% test non liniar transform on bf image and plot
bf_TT = imtransform(bf_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);

figure(1)
imshowpair(bf_T.*2,bf_TT.*2)
hold on
b=A{1};
plot(b(:,2),b(:,1),'-y','LineWidth',2)
plot(input_points(:,1),input_points(:,2),'b*','MarkerSize',20)
plot(base_points(:,1),base_points(:,2),'r*','MarkerSize',20)
legend('Input points','Base points')
title('Transforming to standard fish shape. (Avoid transforming fish shape close to primary tumor)','FontSize',15);
pause(1);
end
end
end

toc; 
disp(' ')

end


% ***INSERTED using CreateNestedCFS_v1_0***
function [ ] = CreateFishStruct_step18_v1_1(specs,aveFish,show)
%% CreateFishStruct_step18
%  Version 1.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 8/24/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%% apply standard fish transform to segmented images
%% Version History:
%  1.1: switch to using dataUpTo and lookup instead of dataAllTimes and config

disp('Step 18: Applying transform functions found in previos step to all images and saving in finals...' );
tic;

% boundary points of standard fish shape
A = bwboundaries(bwmorph(aveFish,'dilate',2));
% get width and heigth of images
ww = specs.maximum_width;
hh = specs.maximum_length;

for fff=BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(fff).lookup;
for t=1:length(lookup)
ff = lookup(t);
% retriavae images from data struct.
body   = BYTIME(t).fishes(ff).rounds(4).body;
spine   = BYTIME(t).fishes(ff).rounds(4).spine;
eye   = BYTIME(t).fishes(ff).rounds(4).eye;
bw   = BYTIME(t).fishes(ff).rounds(4).bw;
bf   = BYTIME(t).fishes(ff).rounds(4).bf;
gfp1   = BYTIME(t).fishes(ff).rounds(4).gfp1;
gfp2   = BYTIME(t).fishes(ff).rounds(4).gfp2;
rfp1   = BYTIME(t).fishes(ff).rounds(4).rfp1;
rfp2   = BYTIME(t).fishes(ff).rounds(4).rfp2;
pure   = BYTIME(t).fishes(ff).rounds(4).pure;
L      = BYTIME(t).fishes(ff).rounds(4).L;
L_divided      = BYTIME(t).fishes(ff).rounds(4).L_divided;

% retriavae transform functions from data struct.
TFORM1 = BYTIME(t).fishes(ff).transforms.lengthScaleTF2;
TFORM2 = BYTIME(t).fishes(ff).transforms.standardFishTF;

% transform images
bf_T    = imtransform(bf,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
bf_TT   = imtransform(bf_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
gfp1_T  = imtransform(gfp1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
gfp1_TT = imtransform(gfp1_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
gfp2_T  = imtransform(gfp2,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
gfp2_TT = imtransform(gfp2_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
rfp1_T  = imtransform(rfp1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
rfp1_TT = imtransform(rfp1_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
rfp2_T  = imtransform(rfp2,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
rfp2_TT = imtransform(rfp2_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
pure_T  = imtransform(pure,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
pure_TT = imtransform(pure_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
body_T  = logical(imtransform(body,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
body_TT = logical(imtransform(body_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
spine_T  = logical(imtransform(spine,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
spine_TT = logical(imtransform(spine_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
eye_T  = logical(imtransform(eye,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
eye_TT = logical(imtransform(eye_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
bw_T  = logical(imtransform(bw,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
bw_TT = logical(imtransform(bw_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
L_T  = uint8(imtransform(L,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
L_TT = uint8(imtransform(L_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
L_divided_T  = uint8(imtransform(L_divided,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
L_divided_TT = uint8(imtransform(L_divided_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));

% save transformed images in finals
BYTIME(t).fishes(ff).finals.body  = body_TT;
BYTIME(t).fishes(ff).finals.spine = spine_TT;
BYTIME(t).fishes(ff).finals.eye   = eye_TT;
BYTIME(t).fishes(ff).finals.bw    = bw_TT;
BYTIME(t).fishes(ff).finals.bf    = bf_TT;
BYTIME(t).fishes(ff).finals.gfp1  = gfp1_TT;
BYTIME(t).fishes(ff).finals.gfp2  = gfp2_TT;
BYTIME(t).fishes(ff).finals.rfp1  = rfp1_TT;
BYTIME(t).fishes(ff).finals.rfp2  = rfp2_TT;
BYTIME(t).fishes(ff).finals.pure  = pure_TT;
BYTIME(t).fishes(ff).finals.L          = L_TT;
BYTIME(t).fishes(ff).finals.L_divided  = L_divided_TT;

% show L_divided_TT in standard fish shape
if show==1
figure(1)
imshow(L_divided_TT,[])
hold on
colormap('jet')
title(['Labeled, clustered and divided segmented image transformed to standard fish shape. Fish #: ' num2str(ff)  ',  Time: ' num2str(t)],'FontSize',13);
b=A{1};
plot(b(:,2),b(:,1),'-y','LineWidth',2)
pause(1);
end
end
end
toc; 
disp(' ');

end


% ***INSERTED using CreateNestedCFS_v1_0***

function [tumors, totals, images] = access_fish_struct_v5_1(name)
%% access_fish_struct_v5_1
%  version 5.1
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 10/11/19
% access_fish_struct_v2_0 takes the struct produced from IdentifyFish2 and
% extracts numerical data characterizing metastasis distribution
%% Version History
%  1.1: not using location anymore as identifier. Just ordering tumors by
%  size assuming the older metastasis are always larger
%  4.0: adding data from intensity image (average intensity of region,
%  intensity stdev, intensity distribution); also adding a solidity measure
%  for the binary image; adding gender id
%  4.1: switch to using dataUpTo and lookup instead of dataAllTimes and config
%  4.2: calculate Absorbance due to tumor from bf
%  5.0: also saving a small subset of images; enabling creation of revised
%  L_divided (L_revised) that has been screened
%  5.1: saving lookup data to csv

tic;
nt_max = length(BYTIME);
nfsh = length(BYTIME(1).dataUpTo{1});

% structure of data in tumors structure array
% tumors(fishID).sizes = [array with sizes of tumors X times]
% tumors(fishID).locs = [array with (xpix,ypix) location of tumor at each time]
% tumors(fishID).reorg = [integer used to represent cluster X times]
% tumors(fishID).count = [number of tumors at each time];
% tumors(fishID).distances(time) = [distance matrix]
% initialize tumors array
tumors = struct('sizes',cell(1,nfsh),'locs',cell(1,nfsh),...
'reorg',cell(1,nfsh),'count',cell(1,nfsh),...
'distances',cell(1,nfsh),'length',cell(1,nfsh),...
'width',cell(1,nfsh),'perimeter',cell(1,nfsh),...
'intensity',cell(1,nfsh),'intensity_std',cell(1,nfsh),...
'convex_solidity',cell(1,nfsh),'gender',cell(1,nfsh),...
'absorbance',cell(1,nfsh),'body_brightness',cell(1,nfsh),...
'absorbance2',cell(1,nfsh),'body_brightness2',cell(1,nfsh));
totals = struct('counts',zeros(1,nt_max),'distrib',zeros(nfsh*4*nt_max,1+nt_max),...
'fish_counts',zeros(1,nt_max),'primary_counts',...
zeros(1,nt_max),'primary_distrib',zeros(nfsh*nt_max,1+nt_max));
images = struct('Ldiv',cell(nt_max,nfsh),'fimg',cell(nt_max,nfsh),...
'bf',cell(nt_max,nfsh),'rgb',cell(nt_max,nfsh),'count',cell(nt_max,nfsh),...
'body',cell(nt_max,nfsh),'spine',cell(nt_max,nfsh),'eye',cell(nt_max,nfsh));
i_distrib = 0;
i_prim = 0;
lookups = zeros(nfsh,nt_max);

% will attempt to extract number and size of metastasis in each fish
for n = BYTIME(1).dataToUse
lookup = BYTIME(1).fishes(n).lookup;
nt = length(lookup);
lookups(n,1:nt) = lookup;
% get final tumor count in fish
tc = 0;
for t=1:nt
nn = lookup(t);
cross_index = BYTIME(t).fishes(nn).cross_index;
tc = tc + size(cross_index,1);
end

% initialize structure array components
tumors(n).sizes = zeros(tc,nt);
tumors(n).locs = zeros(tc,2*nt);
tumors(n).count = zeros(1,nt);
tumors(n).reorg = zeros(tc,nt);
tumors(n).distances = cell(1,nt);
tumors(n).length = zeros(tc,nt);
tumors(n).width = zeros(tc,nt);
tumors(n).perimeter = zeros(tc,nt);
tumors(n).intensity = zeros(tc,nt);
tumors(n).intensity_std = zeros(tc,nt);
tumors(n).convex_solidity = zeros(tc,nt);
tumors(n).gender = BYTIME(1).fishes(n).gender;
tumors(n).forUse = any(BYTIME(1).dataToUse==n);
tumors(n).absorbance = zeros(tc,nt);
tumors(n).body_brightness = zeros(1,nt);
tumors(n).absorbance2 = zeros(tc,nt);
tumors(n).body_brightness2 = zeros(1,nt);

images(1,n).count = nt;
% interate over time points
for i = 1:nt
% want to start from last time point to collect location data
j = nt-i+1;
nn = lookup(j);
% L_divded is an image labeled such that each cluster is a
% different positive integer
Ldiv = BYTIME(j).fishes(nn).rounds(4).L_divided;
% intensity image, scaled assuming int16 image
fimg = double(BYTIME(j).fishes(nn).rounds(4).pure)/65535;
% brightfield image, scaled assuming int16 image
bf = double(BYTIME(j).fishes(nn).rounds(4).bf)/65535;
% green part rgb image, scaled assuming int16 image
rgbg = double(BYTIME(j).fishes(nn).rounds(4).rgb(:,:,2))/65535;
% body and eye bw images
body = BYTIME(j).fishes(nn).rounds(4).body;
eye = BYTIME(j).fishes(nn).rounds(4).eye;
% store number of clusters
tcj = double(max(Ldiv(:)));
tumors(n).count(j) = tcj;
totals.counts(j) = totals.counts(j) + tcj;
%       count the number of fish at each time point
totals.fish_counts(j) = totals.fish_counts(j) + 1;
% count primary tumors (tumors at time t that where present at time 1)
nprim = size(BYTIME(1).fishes(n).cross_index,1);
totals.primary_counts(j) = totals.primary_counts(j) + nprim;

% store images
images(j,n).Ldiv = Ldiv;
images(j,n).fimg = im2uint8(BYTIME(j).fishes(nn).rounds(4).pure);
images(j,n).auto = im2uint8(MergeBF(...
BYTIME(j).fishes(nn).rounds(4).rfp1,BYTIME(j).fishes(nn).rounds(4).rfp2));
%         images(j,n).gfp1 = BYTIME(j).fishes(nn).rounds(4).gfp1;
%         images(j,n).gfp2 = BYTIME(j).fishes(nn).rounds(4).gfp2;
%         images(j,n).rfp1 = BYTIME(j).fishes(nn).rounds(4).rfp1;
%         images(j,n).rfp2 = BYTIME(j).fishes(nn).rounds(4).rfp2;
images(j,n).bf = im2uint8(BYTIME(j).fishes(nn).rounds(4).bf);
images(j,n).rgb = im2uint8(BYTIME(j).fishes(nn).rounds(4).rgb);
images(j,n).body = body;
images(j,n).eye = eye;
images(j,n).spine = BYTIME(j).fishes(nn).rounds(4).spine;

% array to hold the size of each cluster
clus = zeros(tc,1);
% array for locations
locs = zeros(tc,2);
% major axis length
lens = zeros(tc,1);
% minor axis length
wids = zeros(tc,1);
% perimeter
pers = zeros(tc,1);
% intensity
mean_ints = zeros(tc,1);
% intensity_std
int_stds = zeros(tc,1);
% convex_solidity
conv_sol = zeros(tc,1);
% absorbance
absorb = zeros(tc,1);
absorb2 = zeros(tc,1);

% step 1: find cluster sizes and locations using regionprops
stats = regionprops(Ldiv,'Centroid','Area','MajorAxisLength',...
'MinorAxisLength','Perimeter');

% make non-tumor, non-eye region of fish
rest_of_fish = body;
rest_of_fish(bwmorph(logical(Ldiv),'dilate',2)==1) = 0;
rest_of_fish(bwmorph(eye,'dilate',2)==1) = 0;
[~,lastEye] = find(eye==1,1,'last');
rest_of_fish(:,1:lastEye+100) = 0;
other_bf = bf(rest_of_fish==1);
other_rgb = rgbg(rest_of_fish==1);
other_bf_mean = mean(other_bf);
other_rgb_mean = mean(other_rgb);
tumors(n).body_brightness(j) = other_bf_mean;
tumors(n).body_brightness2(j) = other_rgb_mean;
for c = 1:tcj
% get size of cluster
clus(c) = stats(c).Area;
% get location of each cluster
locs(c,:) = stats(c).Centroid;
% dimensions
lens(c) = stats(c).MajorAxisLength;
wids(c) = stats(c).MinorAxisLength;
pers(c) = stats(c).Perimeter;

% analyze intensity image
% 1) identify cluster
reg = return_sub_listL(single(Ldiv),c);
% 2) obtain and analyze intensity data
reg_ints = fimg(reg==c);
if(isempty(reg_ints))
error('Region identification problem at time %i, fish %i',j,nn);
end
mean_ints(c) = mean(reg_ints);
int_stds(c) = std(reg_ints);

% analyze brightfield image for absorbance
reg_bf = bf(reg==c); % brightfield of tumor region
reg_bf_mean = mean(reg_bf);
reg_rgb = rgbg(reg==c);
reg_rgb_mean = mean(reg_rgb);
absorb(c) = log10(other_bf_mean/reg_bf_mean);
absorb2(c) = log10(other_rgb_mean/reg_rgb_mean);

% convex solidity (how much of the convex hull of the region is
% covered by the bw image)
% 1) get convex hull corners
[regy,regx] = find(reg);
% 2) get convex hull bw image
% check colinearity
if(length(regx) > 1)
slopes = (regy(2:end)-regy(1))./(regx(2:end)-regx(1));
isLinear = all(slopes == slopes(1));
else
isLinear = 0;
end
if (~isLinear && length(regx) > 2)
regi = convhull(regx,regy);
reghull = poly2mask(regx(regi),regy(regi),size(Ldiv,1),size(Ldiv,2));
else
regi = 1:length(regx);
reghull = zeros(size(Ldiv));
if(isLinear)
% edge points
if(isinf(slopes(1)))
[~,li] = min(regy);
[~,ri] = max(regy);
x1 = regx(li);
x2 = regx(ri);
y1 = regy(li);
y2 = regy(ri);
else
[~,li] = min(regx);
[~,ri] = max(regx);
x1 = regx(li);
x2 = regx(ri);
y1 = regy(li);
y2 = regy(ri);
end

% draw a line
num_points = 1+round(max(abs(x2-x1),abs(y2-y1)));
xpts = linspace(x1,x2,num_points);
ypts = linspace(y1,y2,num_points);
for ipt = 1:num_points
reghull(ypts(ipt),xpts(ipt)) = 1;
end
end
end
% mark corner points
for ii = 1:length(regi)
reghull(regy(regi(ii)),regx(regi(ii))) = 1;
end
%             reghull = bwmorph(reghull,'thicken',1);
% 3) get solidity as mean value in region covered by convex
% hull
conv_sol(c) = mean(double(Ldiv(reghull==1)))/c;

end

% step 2: reorganize clusters by size;
[~,ord] = sort(clus);
reorg = flipud(ord); % flipped so that largest cluster is first

% step 3: calculate distances between loc and and reference
% calculate distances to all ref locations and store in tcj by tc
% matrix
dist = zeros(tc);
if(j<nt)
ref = tumors(n).locs(:,2*j+1:2*j+2); % loc of clusters at j+1
for c = 1:tcj
for rc = 1:tc
dist(c,rc) = sqrt(sumsqr(locs(c,:)-ref(rc,:)));
end
end
else
ref = locs; % loc of clusters at current time
for c = 1:tcj
for rc = 1:tc
dist(c,rc) = sqrt(sumsqr(locs(c,:)-ref(rc,:)));
end
end
end
% store results
tumors(n).sizes(:,j) = clus;
tumors(n).locs(:,2*j-1:2*j) = locs;
tumors(n).reorg(:,j) = reorg;
tumors(n).distances{j} = dist;
tumors(n).length(:,j) = lens;
tumors(n).width(:,j) = wids;
tumors(n).perimeter(:,j) = pers;
tumors(n).intensity(:,j) = mean_ints;
tumors(n).intensity_std(:,j) = int_stds;
tumors(n).convex_solidity(:,j) = conv_sol;
tumors(n).absorbance(:,j) = absorb;
tumors(n).absorbance2(:,j) = absorb2;

if(tcj > 0)
i_prim = i_prim + 1;
totals.primary_distrib(i_prim,1) = tumors(n).sizes(1,j);
totals.primary_distrib(i_prim,1+j) = 1;
end

for cc = 1:tcj
i_distrib = i_distrib + 1;
totals.distrib(i_distrib,1) = tumors(n).sizes(cc,j);
totals.distrib(i_distrib,1+j) = 1;
end

% step 4: warn about ambiguity when the distances suggest a
% different way to assign clusters than the sizes
[~,mdist] = min(tumors(n).distances{j}(1:tcj,:),[],2);
% mdist = 1:tcj if distance metric agrees with sizes
if(~all(mdist==(1:tcj)'))
disp(['Ambiguity found at fish ' num2str(nn) ' and time ' num2str(j)]);
% display the images to show problem
% 1: load  images
imgj = BYTIME(j).fishes(nn).rounds(4).L_divided;
imgf = BYTIME(nt).fishes(lookup(nt)).rounds(4).L_divided;
%             % 2: reassign interger labels according to reorg
%             reforg = tumors(n).reorg(:,nt);
%             for pix = 1:numel(imgj)
%                 if(imgj(pix) ~= 0)
%                     imgj(pix) = find(reorg==imgj(pix),1);
%                 end
%                 if(imgf(pix) ~= 0)
%                     imgf(pix) = find(reforg==imgf(pix),1);
%                 end
%             end
% 3: mark locations
hwd = 2;
for k = 1:tcj
irj = round(locs(k,1));
icj = round(locs(k,2));
imgj(icj-hwd:icj+hwd,irj-hwd:irj+hwd) = tc+1;
end
for k = 1:tc
irf = round(tumors(n).locs(k,2*nt-1));
icf = round(tumors(n).locs(k,2*nt));
imgf(icf-hwd:icf+hwd,irf-hwd:irf+hwd) = tc+1;
end
% 4: show images
%             figure;
%             subplot(2,1,1)
%             imshow(imgj,[0 tc+1]);
%             subplot(2,1,2)
%             imshow(imgf,[0 tc+1]);
%             colormap('jet');
% %             print(['blank/afs_tumor' num2str(n) '_time' num2str(j) '.jpeg'],'-dpng');
%             pause(1);
% %             close;
% %             disp('done with #4');
end   
end
end

totals.distrib = totals.distrib(1:i_distrib,:);
% some of the primary tumors don't appear at all time points, so I am
% adding these tumors to the distribution at a size smaller than the
% apparent 3 pixel threshold (1 pixel in this case)
true_fish_count = max(totals.fish_counts);
totals.adjusted_primary_count = true_fish_count;
totals.primary_distrib = [1 true_fish_count-totals.fish_counts; totals.primary_distrib(1:i_prim,:)];
totals.distrib = [1 true_fish_count-totals.fish_counts; totals.distrib];
[~,idist_sort] = sort(totals.distrib(:,1));
[~,iprim_sort] = sort(totals.primary_distrib(:,1));
totals.distrib = totals.distrib(idist_sort,:);
totals.primary_distrib = totals.primary_distrib(iprim_sort,:);
% if(totals.distrib(1,1)==0)
%     i_distrib = i_distrib-1;
%     totals.distrib = totals.distrib(2:end,:);
% end

% save results
save([name '_tumors' '.mat'],'tumors');
save([name '_totals' '.mat'],'totals');
save([name '_images' '.mat'],'images');
dlmwrite([name '_lookups.csv'],lookups,'delimiter','\t');
fprintf('Saved files %s.mat\n',name);

toc;

% the inclusion of a division line in L_divided artificially lowers the
% size of two tumors

% should plot clusters over time with labels to demonstrate it worked



end


% ***INSERTED using CreateNestedCFS_v1_0***
%% CreateFishStruct_step19
%  Version 1.2
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 8/24/17
%  Project: Tumor Growth, Logarithmic Continuum Form

function [] = CreateFishStruct_step19_v1_2(name)
%% Version History
%  1.1: save individual fish separately so files are not so large
%  1.2: switch to using dataUpTo and lookup instead of dataAllTimes and config

disp([ 'Step 19: Saving data struct. FISH in current directory as ' name]);
tic;

save_struc(length(BYTIME)).fishes = struct();
for ff = BYTIME(1).dataUpTo{1}
lookup = BYTIME(1).fishes(ff).lookup;
nt = length(lookup);
for t = 1:nt
f = lookup(t);
fields = fieldnames(BYTIME);
for fi = 1:length(fields)
field = fields{fi};
if(strcmp(field,'fishes'))
save_struc(t).fishes = BYTIME(t).fishes(f);
else
save_struc(t).(field) = BYTIME(t).(field);
end
end
end
% save data struct FISH as a .mat file
save_struc = save_struc(1:nt);
save([name '_id' num2str(ff)],'save_struc','-v7.3');
clear save_struc;
save_struc(length(BYTIME)).fishes = struct();
end


toc

end


end