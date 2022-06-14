function [shifted, tdr, tdc] = alignBW_v2_0(img,ref,row,col,max_iter,half_life)
%% alignBW
%  Version 1.3
%  Author: Adeyinka Lesi
%  Date: 9/28/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%  align BW images that have a small shift - not meant to be optimal
%  alignment
%% Version History
%  1.1: img is now can be uncropped
%  1.2: max_iter is an argument
%  1.3: using two stage alignment get better reliability, modified what
%  data types are used - changed alignment algorithm to make it simpler and
%  fix a bug (its so much better that I copied the new alignBW_step to the
%  past versions as well)
%  1.4: half-life put in as a parameter and choosing lowest score candidate
%  as final output
%  2.0: alternative algorithm that looks at 2D regions rather than 1D

% check if either image is blank
if(~any(img(:)) || ~any(ref(:)))
    shifted = img;
    tdr = 0;
    tdc = 0;
    warning('Blank input image');
else
    % caculate initial guess based on centroids of images
%     shift0 = getShiftedImg(img,size(ref),row,col,0,0);
%     [imgys,imgxs] = find(shift0==1);
%     img_cen = [mean(imgxs) mean(imgys)];
%     [refys,refxs] = find(ref==1);
%     ref_cen = [mean(refxs) mean(refys)];
%     dr = round(ref_cen(2)-img_cen(2));
%     dc = round(ref_cen(1)-img_cen(1));
%     [nr,nc] = size(ref);
%     if(dc > 0.1*nc)
%         dc = round(0.1*nc);
%     elseif(dc < -0.1*nc)
%         dc = round(-0.1*nc);
%     end
%     if(dr > 0.1*nr)
%         dr = round(0.1*nr);
%     elseif(dr < -0.1*nr)
%         dr = round(-0.1*nr);
%     end
    dr = 0;
    dc = 0;
    % iterate to improve fit
    MAX_ITER = max_iter;
    score = zeros(1,MAX_ITER+1);
    score_func = @(x) sum(ref(x==1));
    drdc  = zeros(MAX_ITER+1,2);
    score(1) = score_func(getShiftedImg(img,size(ref),row,col,dr,dc));
    drdc(1,:) = [dr dc];
    for n = 1:MAX_ITER
        [ndr, ndc] = alignBW_scored(img,ref,row,col,dr,dc,half_life,score(n),score_func);
        dr = dr + ndr;
        dc = dc + ndc;
        score(n+1) = score_func(getShiftedImg(img,size(ref),row,col,dr,dc));
        drdc(n+1,:) = [dr dc];
        if(ndr^2+ndc^2 == 0)
%                     disp([num2str(n) ' iterations']);
            break;
        elseif(n == MAX_ITER)
%                     disp([num2str(n) ' iterations']);
        end
    end
    [~,imax] = max(score(1:n+1));
    tdr = drdc(imax,1);
    tdc = drdc(imax,2);
    shifted = getShiftedImg(img,size(ref),row,col,tdr,tdc);
%     figure(3);
%     subplot(1,2,1);
%     plot(1:n+1,score(1:n+1));
%     subplot(1,2,2);
%     imshowpair(shifted,ref);
end
end

function [next_dr, next_dc] = alignBW_scored(img,ref,srow,scol,dr,dc,half_life,score,scorer)
    % first stage like before (dr1 and dc1 can be decimals)
    [dr1, dc1] = alignBW_step2D(img,ref,srow,scol,dr,dc,half_life);

    % get step change in both directions
    if(dc1 ~= 0 && dr1 ~= 0)
        if(dc1 > 0)
            idc1 = max(1,round(dc1));
        else
            idc1 = min(-1,round(dc1));
        end
        if(dr1 > 0)
            idr1 = max(1,round(dr1));
        else
            idr1 = min(-1,round(dr1));
        end
    elseif(dc1 ~= 0)
        if(dc1 > 0)
            idc1 = max(1,round(dc1));
        else
            idc1 = min(-1,round(dc1));
        end
        idr1 = 2*round(rand(1))-1;
    elseif(dr1 ~= 0)
        idc1 = 2*round(rand(1))-1;
        if(dr1 > 0)
            idr1 = max(1,round(dr1));
        else
            idr1 = min(-1,round(dr1));
        end
    else
        idc1 = 2*round(rand(1))-1;
        idr1 = 2*round(rand(1))-1;
    end
    
    % score shift
    spdrdc = scorer(getShiftedImg(img,size(ref),srow,scol,dr+idr1,dc+idc1));
    if(spdrdc > score)
        % accept step changes
        next_dr = idr1;
        next_dc = idc1;
    else
        % use interpolation to try and find maximum
        spdr1 = scorer(getShiftedImg(img,size(ref),srow,scol,dr+idr1,dc));
        spdc1 = scorer(getShiftedImg(img,size(ref),srow,scol,dr,dc+idc1));
        spdr2 = scorer(getShiftedImg(img,size(ref),srow,scol,dr-idr1,dc));
        spdc2 = scorer(getShiftedImg(img,size(ref),srow,scol,dr,dc-idc1));
        
        detr = spdr1+spdr2-2*score;
        rel_posr = -0.5*(spdr1-spdr2)/detr;
        if(detr >= 0 || abs(rel_posr) > 1)
            % no local maximum or maximum is outside of range
            if(spdr1 >= spdr2)
                next_dr = idr1;
            else
                next_dr = -idr1;
            end
        else
            % local maximum exists within range
            next_dr = round(idr1*rel_posr);
        end
        detc = spdc1+spdc2-2*score;
        rel_posc = -0.5*(spdc1-spdc2)/detc;
        if(detc >= 0 || abs(rel_posc) > 1)
            % no local maximum or maximum is outside of range
            if(spdc1 >= spdc2)
                next_dc = idc1;
            else
                next_dc = -idc1;
            end
        else
            % local maximum exists within range
            next_dc = round(idc1*rel_posc);
        end
    end
    
end

function [next_dr, next_dc] = alignBW_twoStage(img,ref,srow,scol,dr,dc,half_life)
    % first stage like before (dr1 and dc1 can be decimals)
    [dr1, dc1] = alignBW_step2D(img,ref,srow,scol,dr,dc,half_life);
    % sample function again by changing only one variable
    % figure out new shifts
    if(dc1 ~= 0 && dr1 ~= 0)
        if(dc1 > 0)
            idc1 = max(1,round(dc1));
        else
            idc1 = min(-1,round(dc1));
        end
        if(dr1 > 0)
            idr1 = max(1,round(dr1));
        else
            idr1 = min(-1,round(dr1));
        end
    elseif(dc1 ~= 0)
        if(dc1 > 0)
            idc1 = max(1,round(dc1));
        else
            idc1 = min(-1,round(dc1));
        end
        idr1 = 2*round(rand(1))-1;
    elseif(dr1 ~= 0)
        idc1 = 2*round(rand(1))-1;
        if(dr1 > 0)
            idr1 = max(1,round(dr1));
        else
            idr1 = min(-1,round(dr1));
        end
    else
        idc1 = 2*round(rand(1))-1;
        idr1 = 2*round(rand(1))-1;
    end
    % 2) run alignment
    [dr2x, dc2x] = alignBW_step2D(img,ref,srow,scol,dr,dc+idc1,half_life);
    [dr2y, dc2y] = alignBW_step2D(img,ref,srow,scol,dr+idr1,dc,half_life);
    % 3) get getter steps
    % currently we've made it so that idc1 and idr1 are always non-zero
    % however, if the alignment is insensitive to one both variables, we
    % need to do something different
    difcx = dc2x - dc1;
    difrx = dr2x - dr1;
    difcy = dc2y - dc1;
    difry = dr2y - dr1;
    
    if(difcx==0 && difrx==0 && difcy==0 && difry==0)
        % case 1: all zero
        dx = 0;
        dy = 0;
    elseif(difcx==0 && difrx==0)
        % case 2: insensitive to x
        dx = 0;
        dy = -idr1*(dr1*difry+dc1*difcy)/(difry^2+difcy^2);
    elseif(difcy==0 && difry==0)
        % case 3: insensitive to y
        dx = -idc1*(dc1*difcx+dr1*difrx)/(difcx^2+difrx^2);
        dy = 0;
    elseif(difcx==0 && difcy==0)
        % case 4: The x equation is undefined; need a new equation (here we
        % assume dy is proportional to dx
        dx = -0.5*idc1*dr1/difrx;
        dy = -0.5*idr1*dr1/difry;
    elseif(difrx==0 && difry==0)
        % case 5: The y equation is undefined; need a new equation (here we
        % assume dy is proportional to dx
        dx = -0.5*idc1*dc1/difcx;
        dy = -0.5*idr1*dc1/difcy;
    else
        denom = difcx*difry-difcy*difrx;
        if(denom ~= 0)
            % case 6: both equations are independent
            dx = idc1*(difcy*dr1-difry*dc1)/denom;
            dy = idr1*(difrx*dc1-difcx*dr1)/denom;
        else
            % case 7: equations are ill-defined; need two new equation 
            % (here we assume dy is proportional to dx and we combine the
            % two original equations to make a new one
            dx = -0.25*(dc1*idc1/difcx+dr1*idc1/difrx);
            dy = -0.25*(dc1*idr1/difcy+dr1*idr1/difry);
        end
    end
    
    % prevent too big a step change
    [nr,nc] = size(ref);
    if(dx > 0.25*nc)
        dx = 0.25*nc;
    elseif(dx < -0.25*nc)
        dx = -0.25*nc;
    end
    if(dy > 0.25*nr)
        dy = 0.25*nr;
    elseif(dy < -0.25*nr)
        dy = -0.25*nr;
    end
    
    next_dr = round(dy);
    next_dc = round(dx);
end

function [next_dr, next_dc] = alignBW_step2D(img,ref,srow,scol,dr,dc,half_life)
% distance half life, limits size of next_dr and next_dc (parameter for
% 'body force' used in alignment)

% shift img
shifted = getShiftedImg(img,size(ref),srow,scol,dr,dc);

% get region info
shi_regs = table2array(regionprops('table',shifted,'Area','Centroid'));
ref_regs = table2array(regionprops('table',ref,'Area','Centroid'));
nshi = size(shi_regs,1);
nref = size(ref_regs,1);

% calculate distances between regions in shifted and ref and weights base
% and area and distance
x_dist = zeros(nshi,nref);
y_dist = zeros(nshi,nref);
weights = zeros(nshi,nref);
for i = 1:nshi
    shi_pos = shi_regs(i,2:3);
    shi_area = shi_regs(i,1);
    for j = 1:nref
        ref_pos = ref_regs(j,2:3);
        ref_area = ref_regs(j,1);
        x_dist(i,j) = ref_pos(1)-shi_pos(1);
        y_dist(i,j) = ref_pos(2)-shi_pos(2);
        dist = sqrt(x_dist(i,j)^2+y_dist(i,j)^2);
        weights(i,j) = ref_area*shi_area*exp(-dist*log(2)/half_life);
    end
end

% get weighted average of position shifts
reg_weights = weights./sum(weights,2);
x_pull_reg = sum(x_dist.*reg_weights,2);
y_pull_reg = sum(y_dist.*reg_weights,2);
final_weights = sum(weights.*reg_weights,2); 
final_weights = final_weights/sum(final_weights);
next_dr = final_weights'*y_pull_reg;
next_dc = final_weights'*x_pull_reg;
end

function [next_dr, next_dc] = alignBW_step(img,ref,srow,scol,dr,dc,half_life)
% dimensions (both images should have the same)
[nr,nc] = size(ref);

% distance half life, limits size of next_dr and next_dc (parameter for
% 'body force' used in alignment)

% shift img
shifted = getShiftedImg(img,size(ref),srow,scol,dr,dc);
% get edges for images
edgec_shi = logical([shifted(:,1) shifted(:,2:end)-shifted(:,1:end-1) ...
    -shifted(:,end)]);
edgec_ref = logical([ref(:,1) ref(:,2:end)-ref(:,1:end-1) ...
    -ref(:,end)]);
edger_shi = logical([shifted(1,:); shifted(2:end,:)-shifted(1:end-1,:); ...
    -shifted(end,:)]);
edger_ref = logical([ref(1,:); ref(2:end,:)-ref(1:end-1,:); ...
    -ref(end,:)]);
% variables used for alignment
dr_vals = zeros(1,nc);
dc_vals = zeros(1,nr);
dr_weights = zeros(1,nc);
dc_weights = zeros(1,nr);
% sum_row = sum(shifted,2)'; % used as weights
% sum_col = sum(shifted);

   
% fprintf('Sizes 2: %i,%i,%i,%i\n',size(img,1),size(img,2),...
%     size(shifted,1),size(shifted,2));
% fprintf('Sizes 2: %i,%i,%i,%i\n',size(dif,1),size(dif,2),...
%     size(edgec,1),size(edgec,2));
for row = 1:nr
   % find indices for edges
   nodes_shi = find(edgec_shi(row,:)~=0);
   nodes_ref = find(edgec_ref(row,:)~=0);
   
   [dc_vals(row),dc_weights(row)] = getStepSize(nodes_shi,nodes_ref,half_life);
end

for col = 1:nc
   % find indices for edges
   nodes_shi = find(edger_shi(:,col)~=0);
   nodes_ref = find(edger_ref(:,col)~=0);

   [dr_vals(col),dr_weights(col)] = getStepSize(nodes_shi,nodes_ref,half_life);
end

if(sum(dr_weights) ~= 0)
    dr_weights = dr_weights/sum(dr_weights);
end

if(sum(dr_weights) ~= 0)
    dc_weights = dc_weights/sum(dc_weights);
end
% no need to round here for 2 stage method
next_dr = dr_weights*dr_vals';
next_dc = dc_weights*dc_vals';
% if(isnan(next_dr))
%     pause(1);
% end
end

function [step,weight] = getStepSize(nodes_shi,nodes_ref,half_life)
nnshi = 0.5*length(nodes_shi);
nnref = 0.5*length(nodes_ref);

if(nnshi>0 && nnref>0)
    % edges come in pairs and mark the position of an object
    % will align objects in shifted onto ref by applying a 'body force' to
    % each object in shifted based on the objects in ref. This force is
    % the number of pixels COM of the objects. These forces will be put in a
    % weighted average with weights being an exponential decay term that
    % decreases with distance and depends on the half_life variable - this
    % means the force will be determined by the closest objects
    force = zeros(1,nnshi);
    force_weights = zeros(1,nnshi);
    for nsi = 1:nnshi
        % get position of edges
        ns1 = nodes_shi(2*nsi-1);
        ns2 = nodes_shi(2*nsi);
        % ns_force is the force from each object in ref
        ns_force = zeros(1,nnref);
        % this weight is based on distance from object in shifted as
        % well as size of objects
        ns_weights = zeros(1,nnref);
        for nri = 1:nnref
            % get position of edges
            nr1 = nodes_ref(2*nri-1);
            nr2 = nodes_ref(2*nri);
            % get distance
            dist = 0.5*(nr1-ns1+nr2-ns2);
            ns_force(nri) = dist;
            ns_weights(nri) = (nr2-nr1)*exp(-abs(dist)*log(2)/half_life);
        end
        % force_weights based on size of region and distance
        force_weights(nsi) = (ns2-ns1)*sum(ns_weights);
        ns_weights = ns_weights/sum(ns_weights);
        force(nsi) = ns_force*ns_weights';
    end
    weight = sum(force_weights);
    force_weights = force_weights/sum(force_weights);
    step = force*force_weights';
    weight = weight*exp(-abs(step)*log(2)/half_life);
else
    step = 0;
    weight = 0;
end
end

function [shifted] = getShiftedImg(img,ref_size,row,col,dr,dc)
% if(isnan(dr) || isnan(dc))
%     pause(1);
% end
shifted = false(ref_size);
srow_raw = row-dr;
scol_raw = col-dc;
erow = srow_raw+ref_size(1)-1;
ecol = scol_raw+ref_size(2)-1;
srow = max(1,srow_raw);
scol = max(1,scol_raw);
erow = min(erow,size(img,1));
ecol = min(ecol,size(img,2));

% fprintf('shift_dim: %i, %i; dims: %i,%i,%i,%i\n',shift_row,...
%     shift_col,srow,erow,scol,ecol);
shift_row1 = srow-srow_raw+1;
shift_col1 = scol-scol_raw+1;
shift_row2 = min(ref_size(1),erow-srow+1)+shift_row1-1;
shift_col2 = min(ref_size(2),ecol-scol+1)+shift_col1-1;
shifted(shift_row1:shift_row2,shift_col1:shift_col2) = img(srow:erow,scol:ecol);

end