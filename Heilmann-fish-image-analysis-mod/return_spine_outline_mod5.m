function [final_bw,spine_func] = return_spine_outline_mod5(im_orig,body,body_end,eye_coor,gfp,name)
% This function returns an outline/mask of the fish spine in the tail
% using a bf image of the fish and the body boundary mask
%% Version History
%  mod: changing bwmorph parameters and how endpoint is calculated
%  mod2: it would be really useful to have a more refined backbone outline
%  and there is also the issue that pigmented tumors can interfere with
%  getting an outline
%  mod4: changes to mask_top detection (using body instead of spine_bw) and
%  occulted spine and set limit to amp multiplier for im since when
%  pigmented cells are in frame, value get too high; added name input and
%  outputing image
%  mod5: watch out for spine_bw in two separate pieces

% find last pixel in fish body outline
lastPixel = body_end(1);
% [~,lastPixel] = find(body,1,'last');
        
% we wish to cut out an image of fish tail with spine and muscle (but no background plate)
% at a distance temp_width (150) back from the lastPixel. This cutout will be used to
% find an optimal threshold for distinguishing spine from muscle inside
% fish body
temp_width = 150;

% find upper and lower y values of rectancle for image crop
y_corner     = max([find(body(:,round(lastPixel-temp_width)),1,'first') ,find(body(:,round(lastPixel-temp_width)+100),1,'first')]);
y_corner_low = min([find(body(:,round(lastPixel-temp_width)),1,'last') ,find(body(:,round(lastPixel-temp_width)+100),1,'last')]);

% rectangle for cropping image
spine_rect = [round(lastPixel-temp_width) y_corner+5 +100 y_corner_low-y_corner-10];


amp = 1; % max(1,min(2,65535/crop_limit*.7));
im = im_orig*amp;
crop_bf = imcrop(im,spine_rect);
crop_limit = mean(crop_bf(crop_bf>mean(crop_bf(:))));

% find how dark eye is
bf_eye = double(im(eye_coor(2)+(-10:10),eye_coor(1)+(-10:10)))/65535;
mean_bf_eye = mean(bf_eye(:));
std_bf_eye = std(bf_eye(:));
black_thresh = mean_bf_eye+4*std_bf_eye;

% level for threshold (using Otsu'z method) 
level1 = graythresh(crop_bf(crop_bf<crop_limit));
im_body = double(im(body==1))/65535;
level2 = mean(im_body(im_body>black_thresh));
level = max(level1,level2);
 
spine_bw = 1-imbinarize(im,level); 

% erode body and use this as mask to remove unwanted objects around spine
% outline
body = bwmorph(body,'erode',5);%6
spine_bw(body==0)=0;

% erode spine outline to get rid of small protrusions (blood vessel ect.)
spine_bw = bwmorph(spine_bw,'dilate',2);
spine_bw = bwmorph(spine_bw,'erode',7);

% remove objects smaller than 1000 pixels (noise/dirt)
thres = 1000;
spine_bw = bwareaopen(spine_bw,thres);    

% dilate to smoothen
spine_bw = bwmorph(spine_bw,'dilate',5);
                 
% fill in holes
spine_bw = imfill(spine_bw,'holes');

% smoothen outline
spine_bw = smoothBW(spine_bw,10);  
% make sure there is only one region
labeled = bwlabel(spine_bw);
primary = mode(nonzeros(labeled));
labeled(labeled~=primary) = 0;
spine_bw = logical(labeled);

% spine_bw includes the fish belly and pigmented tumors
% want to use edge detection to find out where the spine of the fish
% actually is throughout the body area

% make mask for location of spine (upper half of fish body)
smoother_spine = smoothBW(spine_bw,35);
end_coor = getBodyEndPoint_v1_0(smoother_spine,20,70);
mask_start = eye_coor(1)+161;
mask_end = end_coor(1);
% modify mask height so there is no noise introduced by the fish belly
% I also need to account for the fact that spine_bw may have some
% protusions that introduce big deviations in to location of the top of the
% spine image
mask_top = zeros(1,mask_end-mask_start+1);
body_wid = zeros(1,mask_end-mask_start+1);
for c = mask_start:mask_end
    r1 = find(body(:,c)==1,1,'first');
    r2 = find(body(:,c)==1,1,'last');
    if(~isempty(r1))
        mask_top(c-mask_start+1) = r1+10;
        body_wid(c-mask_start+1) = r2-r1+1;
    end
end
mask_height = round(max(body_wid)*0.35);

spine_mask = smoother_spine;
spine_mask(end_coor(2)+mask_height:end,:)=0;
spine_mask(:,1:mask_start-1)=0;

% get a smoother body line using a polynomial curve 
imask98 = round(length(mask_top)*0.98); % don't use tail region in fit
nz_top = find(mask_top(1:imask98)~=0);
[topfitx,topfity] = getEvenlySpacedData(nz_top,mask_top(nz_top),20,0.75);
[pfit_top,~,mu_top] = polyfit(topfitx,topfity,3);
smooth_mask_top = round(polyval(pfit_top,1:length(mask_top),[],mu_top));
% remove outliers and get a better fit
top_res = mask_top(nz_top)-smooth_mask_top(nz_top);
std_top_res = std(top_res);
nz_top2 = nz_top(top_res>-std_top_res); % only want to remove stuff that sticks out
[topfitx,topfity] = getEvenlySpacedData(nz_top2,mask_top(nz_top2),20,0.75);
[pfit_top,~,mu_top] = polyfit(topfitx,topfity,4);
smooth_mask_top = round(polyval(pfit_top,1:length(mask_top),[],mu_top));
% calculate mask bottom
mask_bottom = smooth_mask_top + mask_height;
% remove outliers, but not in tail region
max_removal_position = round(eye_coor(1)+0.8*(mask_end-eye_coor(1)));
i_rempos = find(nz_top2+mask_start-1<max_removal_position,1,'last');
std_top = std(mask_top(nz_top2(1:i_rempos)) - ...
    smooth_mask_top(nz_top2(1:i_rempos)));
top_buffer = round(2*std_top);
for c = mask_start:max_removal_position
    spine_mask(1:smooth_mask_top(c-mask_start+1)-1-top_buffer,c) = 0;
end
% remove anything below bottom
for c = mask_start:mask_end
    spine_mask(mask_bottom(c-mask_start+1):end,c) = 0;
end

% use mask to modify image
im_spine = im;
im_spine(spine_mask==0) = 0;

% identify prominent tumors in mask region
gfp_body = double(gfp);
gfp_body(body==0)=0;
[cdf,bin] = histcounts(nonzeros(gfp_body(:)),100,'Normalization','cdf');
gfp_high = min(bin(find(cdf>0.95,1)),7.5e3); % set a maximum to preserve signal for very bright tumors
gfp_body = gfp_body/gfp_high;

% remove dim pixels to create bright_mask
gfp_cutoff = 1;
bright_gfp = gfp_body;
bright_gfp(gfp_body<gfp_cutoff)=0;
bright_mask = bwmorph(logical(bright_gfp),'dilate',10);
% bright_mask = bwmorph(bright_mask,'dilate',11);
bright_mask = imfill(bright_mask,'holes');
gfp_thresh = graythresh(nonzeros(gfp_body));
gfp_bw = imbinarize(gfp_body,gfp_thresh);
% gfp_bw = imbinarize(gfp_body,gfp_thresh*.9);
gfp_bw = bwmorph(gfp_bw,'erode',1);
gfp_bw = bwmorph(gfp_bw,'dilate',11);
gfp_bw = imfill(gfp_bw,'holes');
% gfp_bw = bwmorph(gfp_bw,'erode',5);

% detect really pigmented region (darker than eye)
bf_black = 1-imbinarize(im,black_thresh);
% remove eye area and region behind eye since those are dark and
% will show up
bf_black(:,1:mask_start) = 0;
% clean up binarized image (want to be aggressive with this to
% remove false positives)
bf_black = bwmorph(bf_black,'erode',2);
bf_black = bwareaopen(bf_black,20);
bf_black = bwmorph(bf_black,'dilate',3);
bf_black = imfill(bf_black,'holes');

occulted_spine = gfp_bw;
occulted_spine(bright_mask==0)=0;
occulted_spine(bf_black==0) = 0;
occulted_spine(spine_mask==0)=0;
tooOcculted = false;
if(sum(occulted_spine(:))> 0.5*sum(spine_mask(:)))
    warning('Majority of spine occulted by GFP signal for %s',name);
    tooOcculted = true;
end

% get y gradient
hy = fspecial('sobel');
grad = imfilter(double(im_spine)/double(max(im_spine(:))),hy,'replicate');
% remove edge effect
enum = 2;
spine_mask2 = bwmorph(spine_mask,'erode',enum); 
grad(spine_mask2==0)=0;
% looking for a dark line in the image the middle of that line is zero
% 1) smoothen/blur image to reduce noise
blur_filter = ones(10,3);
grad_blur = imfilter(grad,blur_filter,'replicate');
% 2) find zeros with negative slope line
spine_coor = zeros(1,size(grad_blur,2));
for c = mask_start+enum:mask_end-enum
    r1 = find(spine_mask2(:,c)==1,1,'first');
    r2 = find(spine_mask2(:,c)==1,1,'last');
    line = grad_blur(r1:r2,c);
    sline = sign(line);
    zeros1 = find(sline==0); % value is zero
    zeros2 = find((sline(2:end)-sline(1:end-1))==-2); % transitions from pos to neg
    % put together list of zeros while making sure nothing on the occulted
    % list is included
    occulted = occulted_spine(r1:r2,c);
    zero_list = [zeros1(occulted(zeros1)==0); zeros2(occulted(zeros2)==0)];
    zero_list = sort(zero_list);
    
    % find local slope around each zeros (but skipping zeros close to edge)
    zslope = zeros(length(zero_list),1);
    avg_rad = 10; % averaging radius 1+2*avg numbers used to calculate slope
    for zi = 1:length(zero_list)
        low = max(1,zero_list(zi)-avg_rad);
        high = min(r2-r1+1,zero_list(zi)+avg_rad);
        p = polyfit(low:high,line(low:high)',1);
        zslope(zi) = p(1);
    end
    % select point
    % want maximum negative slope
    slope_mag = max(0,-zslope);
    % only search look at the first few zeros
    limit = min(5,length(slope_mag));
    if(~isempty(zero_list))
        [~,i_max] = max(slope_mag(1:limit));
        spine_coor(c) = zero_list(i_max)+r1-1;
    end
end

nz_coor = find(spine_coor~=0);
badSpineDetection = false;
if(isempty(nz_coor))
    error('No spine points detected for %s',name);
elseif(length(nz_coor) < 0.5*length(mask_top))
    warning('Less than half of the spine points were detected for %s',name);
    badSpineDetection = true;
end
% spine_coor_avg = zeros(1,length(nz_coor));
nzi_lim = 5;
dis_lim = 10;
spine_coor_avg = getSelectRunningAverage(spine_coor,nz_coor,nzi_lim,dis_lim);
% for i = 1:length(nz_coor)
%     low1 = max(1,i-nzi_lim);
%     low2 = find(nz_coor(i)-nz_coor(low1:i)<dis_lim,1)+low1-1;
%     high1 = min(length(nz_coor),i+5);
%     high2 = find(nz_coor(i:high1)-nz_coor(i)<dis_lim,1,'last')+i-1;
%     spine_coor_avg(i)=mean(spine_coor(nz_coor(low2:high2)));
% end

% polynomial fit to remove outliers
% should use a subset of points for fit
nfit_max = 100;
% nfit = min(nfit_max,nz_coor(end)-nz_coor(1)+1);
% spacing = (nz_coor(end)-nz_coor(1)+1)/(nfit-1);
% min_spacing = 0.75*spacing;
[fitx,fity] = getEvenlySpacedData(nz_coor,spine_coor_avg,nfit_max,0.75);

% fitx = zeros(1,nfit);
% fity = zeros(1,nfit);
% last_used = 1;
% fitx(1) = nz_coor(last_used);
% fity(1) = spine_coor_avg(last_used);
% nfit_final = 1;
% for i = 2:nfit
%     high1 = find(nz_coor(last_used+1:end)-nz_coor(1)+1>=spacing*(i-1),1);
%     high2 = find(nz_coor(last_used+1:end)-nz_coor(last_used)>=min_spacing,1);
%     high = max(high1,high2);
%     if(isempty(high))
%         % ran out of points
%         break;
%     elseif(high==1)
%         last_used = last_used+1;
%     else
%         if(spacing*(i-1)+1>0.5*(nz_coor(high)+nz_coor(high-1)))
%             last_used = last_used+high;
%         else
%             last_used = last_used+high-1;
%         end
%     end
%     fitx(i) = nz_coor(last_used);
%     fity(i) = spine_coor_avg(last_used);
%     nfit_final = nfit_final + 1;
% end
% fitx = fitx(1:nfit_final);
% fity = fity(1:nfit_final);      

[pfit_coef,~,mu1] = polyfit([eye_coor(1) fitx body_end(1)],...
    [eye_coor(2) fity body_end(2)],3);
spine_fit = polyval(pfit_coef,1:length(spine_coor),[],mu1);

% get residuals, and standard dev of residuals
res = spine_fit(nz_coor)-spine_coor(nz_coor);
std_res = std(res);
nz_coor2 = nz_coor(abs(res)<2*std_res);

% get better spine line
spine_coor_avg2 = getSelectRunningAverage(spine_coor,nz_coor2,nzi_lim,dis_lim);
nfit_max2 = 100;
[fit2x,fit2y] = getEvenlySpacedData(nz_coor2,spine_coor_avg2,nfit_max2,0.75);

[pfit_coef2,S2,mu2] = polyfit([eye_coor(1) fit2x body_end(1)],...
    [eye_coor(2) fit2y body_end(2)],4);
[spine_fit2,sfdel2] = polyval(pfit_coef2,1:length(spine_coor),S2,mu2);

spine_func = struct('observed_x',nz_coor,'observed',spine_coor_avg2,...
    'coefs',pfit_coef2,'scaler',mu2,'error',sfdel2,'calculated_x',...
    eye_coor(1):body_end(1),'calculated',spine_fit2(eye_coor(1):body_end(1)));

final_bw = spine_bw;
spine_rad = 9;
% paint a stripe centered along the spine
for c = nz_coor2(1):nz_coor2(end)
    spp = round(spine_fit2(c));
    final_bw(spp-spine_rad:spp+spine_rad,c) = 1;
    final_bw(1:spp-spine_rad-1,c) = 0;
    final_bw(spp+spine_rad+1:end,c) = 0;
end
max_col = nz_coor2(end);
final_bw(:,max_col+1:end)=0;

% figure(2);
% imshowpair(imresize(spine_bw+bf_black,0.5),imresize(final_bw,0.5));

fig3 = figure('Visible','off');
subplot(2,2,1);
imshow(im);
title(['Amp = ' num2str(amp) ', Thresh = ' num2str(level)]);
hold on;
plot(nz_coor,spine_coor_avg,'LineWidth',1);
plot(nz_coor2,spine_coor_avg2,'LineWidth',2);
plot(mask_start:mask_end,[smooth_mask_top-top_buffer;mask_bottom]);
hold off;
subplot(2,2,2);
imshow(abs(grad_blur)+double(occulted_spine)*0.25*max(abs(grad_blur(:))),[]);
if(tooOcculted)
    title('Too occulted');
end
hold on;
% plot(nz_coor,spine_coor(nz_coor));
% plot(nz_coor,spine_coor_avg);
% plot(1:length(spine_coor),spine_fit,'g-',1:length(spine_coor),spine_fit2,'r-');
plot(1:length(spine_coor),spine_fit2,'b-');
plot(eye_coor(1),eye_coor(2),'o',body_end(1),body_end(2),'o');
hold off;
subplot(2,2,3);
plot(nz_coor,spine_coor(nz_coor),nz_coor,spine_coor_avg,1:length(spine_coor),spine_fit,fitx,fity,'o');
axis([mask_start mask_end min(fity)-10 max(fity)+10]);
if(badSpineDetection)
    title('Bad Detection');
end
subplot(2,2,4);
plot(nz_coor2,spine_coor(nz_coor2),nz_coor2,spine_coor_avg2,1:length(spine_coor),spine_fit2,fit2x,fit2y,'o');
axis([mask_start mask_end min(fity)-10 max(fity)+10]);
slash_loc = find((name=='/')+(name=='\'),1,'last');
if(isempty(slash_loc))
    slash_loc = 0;
end
% stv = suptitle(name(slash_loc+1:end));
% stv.Interpreter = 'none';
title(name(slash_loc+1:end),'interpreter','none');
saveas(fig3,name,'jpg');
close;
% pause(.1);
end

function [select_avg] = getSelectRunningAverage(data,selected,ind_lim,dis_lim)
% calculate a running average - but the data is variably spaced apart based
% on the selected array of indices
select_avg = zeros(1,length(selected));
for i = 1:length(selected)
    low1 = max(1,i-ind_lim);
    low2 = find(selected(i)-selected(low1:i)<dis_lim,1)+low1-1;
    high1 = min(length(selected),i+5);
    high2 = find(selected(i:high1)-selected(i)<dis_lim,1,'last')+i-1;
    select_avg(i)=mean(data(selected(low2:high2)));
end
end

function [fitx,fity] = getEvenlySpacedData(xdat,ydat,nfit_max,min_spacing_frac)
% select at maximum nfit_max points out of data to be used for polynomial
% curve fitting

nfit = min(nfit_max,xdat(end)-xdat(1)+1);
spacing = (xdat(end)-xdat(1)+1)/(nfit-1);
min_spacing = min_spacing_frac*spacing;

fitx = zeros(1,nfit);
fity = zeros(1,nfit);
last_used = 1;
fitx(1) = xdat(last_used);
fity(1) = ydat(last_used);
nfit_final = 1;
for i = 2:nfit
    high1 = find(xdat(last_used+1:end)-xdat(1)+1>=spacing*(i-1),1);
    high2 = find(xdat(last_used+1:end)-xdat(last_used)>=min_spacing,1);
    high = max(high1,high2);
    if(isempty(high))
        % ran out of points
        break;
    elseif(high==1)
        last_used = last_used+1;
    else
        if(spacing*(i-1)+1>0.5*(xdat(high)+xdat(high-1)))
            last_used = last_used+high;
        else
            last_used = last_used+high-1;
        end
    end
    fitx(i) = xdat(last_used);
    fity(i) = ydat(last_used);
    nfit_final = nfit_final + 1;
end
fitx = fitx(1:nfit_final);
fity = fity(1:nfit_final); 
end