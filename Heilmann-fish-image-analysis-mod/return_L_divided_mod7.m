function [L_divideds,V_merges,cross_index,cross_links]  = return_L_divided_mod7(objs)
% returnL_divided returns three labeled matrix where metastases and primary tumor have been seperated (if the had merged previos to time t)
%% Version History
%  simplified code to prevent an error (the loop that checks if the
%  clusters in rClist are present in rLlist and adds a dot if not)
%  8/21/17: modified bwmask calculation to draw a line when points are
%  colinear
%  mod6 fixes squiggly lines from voronoi method separation
%  mod7 changes the Bad association error to a warning, removes bwmask,
%  adds extra Voronoi points for virgin tumor area so that area gets
%  segmented more reasonably
%  2/15/20: fixed bug caused by holes in regions (used imfill on L_k and
%  L_k_combo)

nt = length(objs);
L_divideds = cell(1,nt);
V_merges = cell(1,nt-1);
V_merge=[];

[~,index] = return_met_events_mod3(objs);
% the number corresponding to the same cluster at the next time frame
% cross_links{1}(2) is the number of cluster 2 at time 1 at time 2
cross_links = cell(1,nt-1);
% list of the number corresponding to each new metastasis at each L_divided
% image the numbers for each metastasis start at the time the metastasis
% first appears
cross_index = cell(1,nt); 

% iterate through objs to get BW image and distance matrix
nmets = zeros(1,nt);
for t = 1:nt
    objs(t).BW = logical(objs(t).L_t1);
    objs(t).L_labeled = bwlabel(objs(t).L_t1);
    objs(t).DM = return_dist_matrix(objs(t).BW);
    if(t>1)
        nmets(t) = nmets(t-1) + length(index{t});
    else
        nmets(t) = length(index{t});
    end
    cross_index{t} = zeros(length(index{t}),nt-t+1);
end

% a dilation number used in the code
DC = zeros(1,nt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D1
% L_D1 never needs to be seperated...
L_divideds{1} = objs(1).L_t1;
cross_index{1}(:,1) = index{1}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D7
L_curDivided = 0;
for t = 2:nt   
    disp(objs(t).fishID);
    % find out it L_D7 need to be seperated
    %     cmlist = zeros(nmets(t),2);
    %     cpos = 1;
    %     for tt = 1:t
    %         cmlist(cpos:cpos+length(index{tt}),:) = cmmets{tt};
    %         cpos = cpos+length(index{tt})+1;
    %     end
    L_prevDivided = L_curDivided;
    
%     tempL1 = L_divideds{t-1};
%     tempL2 = objs(t).L_t1;
%     tempL2(tempL1==0)=0; % review any new clusters
    
%     % use centroids to match clusters accross time
%     rL = getCentroids_v1_1(objs(t).L_t1);   % get CM's of blobs in this time point
%     rC = getCentroids_v1_1(L_divideds{t-1}); % get CM's of objects in previous time point

%     % list of indexes of blops in labeled image
%     listBlobsNotMerge = zeros(1,length(rC)); % length(rC) is number of clusters/objects at t-1
% 
%     % for each object at t-1 find the closest blob at t
%     for k=1:length(rC)
%         dist = 1e10;
%         index1 = 0;
%         for n=1:length(rL)
%             
%             cmL = rL(n).Centroid; % (Note 'Centroid' returns first x = col then y = row)
%             cmC = rC(k).Centroid; % (Note 'Centroid' returns first x = col then y = row)
%             
%             temp_dist = sqrt((cmL(1)-cmC(1))^2 + (cmL(2)-cmC(2))^2 ); % distance between n'th blop at time+t and k'th object at time=t-1
%             
%             if dist>temp_dist && temp_dist<55 % save smallest distance
%                 dist = temp_dist;
%                 index1 = n;
%             end
%         end
%         listBlobsNotMerge(k) = index1; % put index of closest blop D7 in list
%     end
        

   

    % if number met events (S(1)) is 1 or if number of indexses is the same for
    % L_D1 and L_D7(L_D1==0)=0 then L_D7 needs no division
    if  nmets(t)<=1 
        L_curDivided = 0;
        if(max(objs(t).L_t1(:))>1)
            warning(['Separate clusters grouped as one tumor with ' objs(t).fishID]);
        end
        L_divideds{t} = double(logical(objs(t).L_t1)); % make sure this is clustered....
    
%     elseif  length(listBlobsNotMerge)==length(unique(listBlobsNotMerge)) && isempty(find(listBlobsNotMerge==0,1))==1 % if number of unique blop indexes is same as number of object no blop has merged!
%     
%         [L_divideds{t},numMetEvents ] = return_clustered_C(objs(t).BW,return_dist_matrix(objs(t).BW),20,listBlobsNotMerge);
%         L_curDivided = 1;
        
    else
        L_curDivided = 1;
        objects = struct('voronoiPointsX',cell(1,nmets(t)),'voronoiPointsY',cell(1,nmets(t)),...
                         'L_combo',cell(1,nmets(t)),'L_bounds',cell(1,nmets(t)),...
                         'L_props',cell(1,nmets(t)));
        objectCounter = 0;
        pointCounter = 0;
        
        [rows,cols] = size(objs(t).L_t1);
        
        % iterate through metastatic events which happened before
        
        % make a consolidated image for each distinct tumor. make sure
        % there is no overlap
        L_covered = zeros(size(objs(t).L_t1));
        for tt = 1:t-1
            for k=1:length(index{tt})
                objectCounter=objectCounter+1;
                L_k_combo = false(size(objs(t).L_t1));
                for ii = tt:t-1
                    L_k = return_sub_listL(L_divideds{ii},cross_index{tt}(k,ii-tt+1));
                    L_k_combo(L_k~=0) = 1;
                end
                L_k_combo = imfill(L_k_combo,'holes');
                objects(objectCounter).L_combo = L_k_combo;
                objects(objectCounter).L_bounds = bwboundaries(L_k_combo);
                objects(objectCounter).L_props = regionprops(L_k_combo,'Area');
                L_covered = L_covered+L_k_combo*2^(objectCounter-1);
            end
        end
        
        % make images of new metastasis
        for k=1:length(index{t})
            objectCounter=objectCounter+1;
            L_k = logical(return_sub_listL(objs(t).L_t1,index{t}(k)));
            if L_prevDivided==0 % not sure this is necessary anymore
                L_k = bwmorph(L_k,'dilate',DC(t)); 
            end
            L_k = imfill(L_k,'holes');
            objects(objectCounter).L_combo = L_k;
            objects(objectCounter).L_bounds = bwboundaries(L_k);
            objects(objectCounter).L_props = regionprops(L_k,'Area');
            % there shouldn't be any overlap with any other regions
            L_covered = L_covered+L_k*2^(objectCounter-1);
        end
        
        % rectify L_combo images so there is no ambiguity between regions
        % L_covered is labeled such that any combination of overlaps is
        % distinguishable using binary math
        base_lab = 2.^(0:objectCounter-1);
        lc_overlaps = setdiff(nonzeros(unique(L_covered)),base_lab);
        overlap_comp = cell(1,length(lc_overlaps));
        for ll = 1:length(lc_overlaps)
            used_lab = false(1,length(base_lab));
            compsum = lc_overlaps(ll);
            for ii = 1:length(base_lab)
                iii = length(base_lab)-ii+1;
                if(compsum>=base_lab(iii))
                    used_lab(iii) = true;
                    compsum = compsum-base_lab(iii);
                end
            end
            overlap_comp{ll} = find(used_lab);
        end
        
        % get boundary of overlap region and use to divide up the overlap
        % region appropriately
        lc_areas = zeros(1,objectCounter);
        lc_bounds = cell(1,objectCounter);
        for k = 1:objectCounter
            for kk = 1:length(objects(k).L_props)
                lc_areas(k) = lc_areas(k)+objects(k).L_props(kk).Area;
            end
            B1 = objects(k).L_bounds;
            nB1 = 0;
            for ib = 1:length(B1)
                nB1 = nB1+length(B1{ib});
            end
            B1lin = zeros(1,nB1);
            nb = 0;
            for ib = 1:length(B1)
                B1i = B1{ib};
                B1lin(nb+(1:size(B1i,1))) = rows*(B1i(:,2)-1)+B1i(:,1);
                nb = nb+size(B1i,1);
            end
            lc_bounds{k} = B1lin;
        end
        
        remLc = false([size(objs(t).L_t1) objectCounter]);
        for ll = 1:length(lc_overlaps)
            Lover = logical(return_sub_listL(L_covered,lc_overlaps(ll)));
            Bover = bwboundaries(Lover);
            if(length(overlap_comp{ll})==2)
                olr1 = overlap_comp{ll}(1);
                olr2 = overlap_comp{ll}(2);
                % olr_out = [];
            else
                %(length(overlap_comp{ll})>2)
                % for simplicity will divide region between two biggest
                % regions
                ol_areas = lc_areas(overlap_comp{ll});
                [~,ol_ord] = sort(ol_areas,'descend');
                olr1 = overlap_comp{ll}(ol_ord(1));
                olr2 = overlap_comp{ll}(ol_ord(2));
                % olr_out = setdiff(overlap_comp{ll},[olr1, olr2]);
            end
            B1 = lc_bounds{olr1};
            B2 = lc_bounds{olr2};
            
            % 1) divide boundary of overlap regions by which Lc boundary
            % region they are part of
            % 2) find a concensus boundary using Lc properties as a factor
            % 3) modify Lc regions with new info
            for bo = 1:length(Bover)
                olB = Bover{bo};
                nolB = size(olB,1);
                Blin = rows*(olB(:,2)-1)+olB(:,1);
                iol1 = find(ismember(Blin,B1));
                iol2 = find(ismember(Blin,B2));
                % assign start and end points of boundary region
                % start will be a point in both B1 and B2
                % end will either be a second point in both B1 and B2 (case
                % of two overlapping regions) or some point with no overlap
                % with B1 and B2 (case of multiple overlapping regions
                % where those other regions have been ignored)
                inter_iol = intersect(iol1,iol2);
                diff_iol = setdiff(1:nolB-1,[iol1;iol2]); % last olB point is a duplicate
                if(length(inter_iol)==2)
                    starti = inter_iol(1);
                    endi = inter_iol(2);
                elseif(length(inter_iol)>2)
                    % pick points furthest apart
                    starti = inter_iol(1);
                    disti = min(abs(inter_iol(2:end)-starti),...
                        nolB-1-abs(inter_iol(2:end)-starti));
                    [~,idi] = max(disti);
                    endi = inter_iol(idi+1);
                elseif(length(inter_iol)==1)
                    % pick end point as furthest point from start that does
                    % not overlap either region, or otherwise, the furthest
                    % point in general
                    starti = inter_iol(1);
                    if(~isempty(diff_iol))
                        end_cands = diff_iol;
                    else
                        end_cands = 1:nolB-1;
                    end
                    disti = min(abs(end_cands-starti),...
                        nolB-1-abs(end_cands-starti));
                    [~,idi] = max(disti);
                    endi = end_cands(idi);
                else
                    % no intersects
                    if(~isempty(iol1)&&~isempty(iol2))
                        % pick based on the endpoints of iol1 and iol2
                        starti = iol1(1);
                        end_cands = iol2;
                        disti = min(abs(end_cands-starti),...
                            nolB-1-abs(end_cands-starti));
                        [~,idi] = max(disti);
                        endi = end_cands(idi);
                    elseif(~isempty(iol1))
                        % pick solely based on iol1
                        starti = iol1(1);
                        end_cands = iol1;
                        disti = min(abs(end_cands-starti),...
                            nolB-1-abs(end_cands-starti));
                        [~,idi] = max(disti);
                        endi = end_cands(idi);
                    else
                        % pick solely based on iol2
                        starti = iol2(1);
                        end_cands = iol2;
                        disti = min(abs(end_cands-starti),...
                            nolB-1-abs(end_cands-starti));
                        [~,idi] = max(disti);
                        endi = end_cands(idi);
                    end
                end
                % one curve goes clockwise from starti to endi while the
                % other goes counter clockwise; however I want to the same
                % number of points in both curves so one will have
                % duplicates
                if(starti>endi) % easier to have known direction
                    tempi = starti;
                    starti = endi;
                    endi = tempi;
                end
                maxlen = max(abs(endi-starti)+1,nolB-abs(endi-starti));
                cwi = round(linspace(starti,endi,maxlen));
                ccwi = round(linspace(nolB-1+starti,endi,maxlen));
                ccwi(ccwi>=nolB) = ccwi(ccwi>=nolB)-nolB+1;
                % get weighted average row and col positions from start to end
                cw_B1_count = sum(ismember(cwi,iol1));
                cw_B2_count = sum(ismember(cwi,iol2));
                isCW1 = cw_B1_count>=cw_B2_count; % indicates if cw curve is part of region 1's interface
                if(isCW1)
                    intwei = lc_areas([olr1 olr2]);
                else
                    intwei = lc_areas([olr2 olr1]);
                end
                intwei = intwei/sum(intwei);
                cwr = olB(cwi,1);
                cwc = olB(cwi,2);
                ccwr = olB(ccwi,1);
                ccwc = olB(ccwi,2);
                mr = round([cwr ccwr]*intwei');
                mc = round([cwc ccwc]*intwei');
                
                % imshow(L_covered,[]);
                % hold on; plot(cwc,cwr,'*',ccwc,ccwr,'*',mc,mr,'o');
                
                % mark areas that need to be removed
                remL1 = false(size(objs(t).L_t1));
                remL2 = false(size(objs(t).L_t1));
                if(isCW1)
                    for ii = 1:maxlen
                        remL1(mr(ii),mc(ii)) = true;
                        remL2(mr(ii),mc(ii)) = true;
                        remL1(cwr(ii),cwc(ii)) = true;
                        remL2(ccwr(ii),ccwc(ii)) = true;
                    end
                else
                    for ii = 1:maxlen
                        remL1(mr(ii),mc(ii)) = true;
                        remL2(mr(ii),mc(ii)) = true;
                        remL2(cwr(ii),cwc(ii)) = true;
                        remL1(ccwr(ii),ccwc(ii)) = true;
                    end
                end
                % close areas outlined by boundary points
                remL1 = imfill(bwmorph(remL1,'close'),'holes');
                remL2 = imfill(bwmorph(remL2,'close'),'holes');
                remLc(:,:,olr1) = remLc(:,:,olr1)|remL1;
                remLc(:,:,olr2) = remLc(:,:,olr2)|remL2;
            end
        end
        
        % remodel L_combos
        L_covered2 = zeros(size(objs(t).L_t1));
        L_overlap = zeros(size(objs(t).L_t1));
        interface_threshold = floor(objs(t).threshold/2);
        for k = 1:objectCounter
            objects(k).L_combo = objects(k).L_combo&(~remLc(:,:,k));
            objects(k).L_bounds = bwboundaries(objects(k).L_combo);
            objects(k).L_props = regionprops(objects(k).L_combo,'Area');
            L_covered2 = L_covered2+objects(k).L_combo*2^(k-1);
            L_overlap = L_overlap+bwmorph(objects(k).L_combo,'dilate',interface_threshold);
        end
        
        % identify regions of overlap
        L_interface = L_overlap>1;
        
        % get boundary voronoi points
        [bdxx,bdyy] = return_boundary_voronoi_points(L_interface,interface_threshold);
        bdlin = sub2ind(size(L_interface),bdyy,bdxx);
        
        %         if(~isempty(lc_overlaps))
        %             imshow(L_covered2+L_interface,hot(max(max(L_covered2+L_interface))));
        %             hold on;
        %             plot(bdxx,bdyy,'o');
        %             pause(1);
        %         end
        
        % get voids around clusters
        allmask = false(size(objs(t).BW));
        for ttt = 1:t
            allmask(objs(ttt).BW==1)=1;
        end
        
        % get Voronoi points from saved images
        % Voronoi points for different objects may be too close - this may
        % cause problems down the line with segmentation. I am making a new
        % voronoi point get-function that takes other objects into account
        
        for k = 1:objectCounter
            % find boundary points that lie on L_combo
            [vxx,vyy] = return_object_voronoi_points_mod5(objects(k).L_combo.*allmask,objs(t).threshold,objects(1:k-1));
            bdonLc = find(objects(k).L_combo(bdlin));
            toAdd = true(1,length(bdonLc));
            for ii = 1:length(bdonLc)
                bdi = bdonLc(ii);
                for jj = 1:length(vxx)
                    dist = sqrt((vxx(jj)-bdxx(bdi))^2+(vyy(jj)-bdyy(bdi))^2);
                    if dist<interface_threshold
                        toAdd(ii) = false;
                    end
                end
            end
            objects(k).voronoiPointsX = [vxx;bdxx(bdonLc(toAdd))];
            objects(k).voronoiPointsY = [vyy;bdyy(bdonLc(toAdd))];
            pointCounter = pointCounter+length(vxx)+sum(toAdd);
        end
        
        
        %         for tt = 1:t-1
        %             for k=1:length(index{tt})
        %                 % the time points between tt and t can also provide
        %                 % information on the cluster, so points from that should
        %                 % also be found
        %                 L_k_combo = false(size(objs(t).L_t1));
        %                 for ii = tt:t-1
        %                     L_k = return_sub_listL(L_divideds{ii},cross_index{tt}(k,ii-tt+1));
        %                     L_k = bwmorph(L_k,'dilate',DC(tt));
        %                     L_k_combo(L_k~=0) = 1;
        %                 end
        %                 [vxx,vyy] = return_object_voronoi_points_mod5(L_k_combo.*allmask,objs(t).threshold,objects(1:objectCounter-1));
        %                 imshow(L_k_combo);
        %                 hold on;
        %                 plot(vxx,vyy,'*');
        %                 hold off;
        %                 pause(1);
        %
        %                 % weed out points which are too close
        %                 % if points are too close replace one of them with (0,0) - then later remove
        %                 for ii=1:length(vxx)
        %                     if(vxx(ii) > 0)
        %                         for jj=ii+1:length(vxx)
        %                             dist = sqrt((vxx(ii)-vxx(jj))^2 +(vyy(ii)-vyy(jj))^2);
        %                             if dist < objs(t).threshold
        %                                 vxx(jj)= 0;
        %                                 vyy(jj)= 0;
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %                 % remove zeros!
        %                 vxx = round(nonzeros(vxx));
        %                 vyy = round(nonzeros(vyy));
        %
        %                 objects(objectCounter).voronoiPointsX = vxx;
        %                 objects(objectCounter).voronoiPointsY = vyy;
        %                 pointCounter = pointCounter+length(vxx);
        %             end
        %         end
        
        %         % iterate through new metastasis
        %         for k=1:length(index{t})
        %             objectCounter=objectCounter+1;
        %
        %             L_k = return_sub_listL(objs(t).L_t1,index{t}(k));
        %             if L_prevDivided==0
        %                 L_k = bwmorph(L_k,'dilate',DC(t));
        %             end
        %
        %             [xx,yy]= return_object_voronoi_points_mod5(L_k,objs(t).threshold,objects(1:objectCounter-1));
        %
        %             %%%%
        %             % weed out points which are to close
        %             % if points are too close replace one of them with (0,0) - then later remove
        %             for ii=1:length(xx)
        %                 if(xx(ii) > 0)
        %                     for jj=ii+1:length(xx)
        %                         dist = sqrt((xx(ii)-xx(jj))^2 +(yy(ii)-yy(jj))^2);
        %                         if dist < objs(t).threshold
        %                             xx(jj)= 0;
        %                             yy(jj)= 0;
        %                         end
        %                     end
        %                 end
        %             end
        %
        %             % remove zeros!
        %             xx = round(nonzeros(xx));
        %             yy = round(nonzeros(yy));
        %
        %             objects(objectCounter).voronoiPointsX = xx;
        %             objects(objectCounter).voronoiPointsY = yy;
        %             pointCounter = pointCounter + length(xx);
        %         end
        
        % add in Voronoi points for new area in old tumors - during the
        % segmentation region later, the areas will be parcelled out in
        % smaller chunks and the results will look more reasonable
        newGrowth = objs(t).BW;
        oldBW = bwmorph(logical(L_covered),'dilate',objs(t).threshold);
        newGrowth(oldBW) = 0;
        [ngxx,ngyy]= return_object_voronoi_points_mod5(newGrowth,objs(t).threshold,objects(1:objectCounter));
        
        % extract all points into one vector so they can be feed to next function
        
        x = zeros(pointCounter+length(ngxx),1);
        y = zeros(pointCounter+length(ngyy),1);
        npt = 0;
        for k=1:length(objects)
            npt_k = length(objects(k).voronoiPointsX);
            x(npt+(1:npt_k)) = objects(k).voronoiPointsX;
            y(npt+(1:npt_k)) = objects(k).voronoiPointsY;
            npt = npt+npt_k;
        end
        x(pointCounter+1:end) = ngxx;
        y(pointCounter+1:end) = ngyy;
        
        %%%%
        
        %         % create a convex hull mask that covers all visible clusters
        %         [bwy1,bwx1] = find(objs(t).BW);
        %         bwx = [bwx1;x];
        %         bwy = [bwy1;y];
        %         % check colinearity
        %         if(length(bwx) > 1)
        %             slopes = (bwy(2:end)-bwy(1))./(bwx(2:end)-bwx(1));
        %             isLinear = all(slopes == slopes(1));
        %         else
        %             isLinear = 0;
        %         end
        %         if(~isLinear && length(bwx)> 2)
        %             % only use convhull if there are 3 or more vertices
        %             bwi = convhull(bwx,bwy);
        %             bwmask = poly2mask(bwx(bwi),bwy(bwi),rows,cols);
        %         else
        %             bwi = 1:length(bwx);
        %             bwmask = zeros(rows,cols);
        %             if(isLinear)
        %                    % edge points
        %                    if(isinf(slopes(1)))
        %                        [~,li] = min(bwy);
        %                        [~,ri] = max(bwy);
        %                        x1 = bwx(li);
        %                        x2 = bwx(ri);
        %                        y1 = bwy(li);
        %                        y2 = bwy(ri);
        %                    else
        %                        [~,li] = min(bwx);
        %                        [~,ri] = max(bwx);
        %                        x1 = bwx(li);
        %                        x2 = bwx(ri);
        %                        y1 = bwy(li);
        %                        y2 = bwy(ri);
        %                    end
        %                    % draw a line
        %                    num_points = 1+round(max(abs(x2-x1),abs(y2-y1)));
        %                    xpts = round(linspace(x1,x2,num_points));
        %                    ypts = round(linspace(y1,y2,num_points));
        %                    for ipt = 1:num_points
        %                        bwmask(ypts(ipt),xpts(ipt)) = 1;
        %                    end
        %             end
        %         end
        %         % fill in key points (in case poly2mask misses some things)
        %         for bwii = 1:length(bwi)
        %             bwmask(bwy(bwi(bwii)),bwx(bwi(bwii))) = 1;
        %         end
        %         bwmask = bwmorph(bwmask,'dilate',objs(t).threshold);
        
        % get segmented regions based on voronoi algorithm
        V_labeled = return_labeled_voronoi_mod2(x,y,rows,cols,objs(t).threshold);
        
        % use convex hull mask to confine V_labeled to regions were objects are
        % (as is, don't need bwmask; may be deprecated later)
        %         V_labeled(bwmask==0) = 0;
        % the dark mask adds the separations in real image to V_labeled.
        % The purpose of this is to help correctly associate regions with
        % the appropriate objects. Using Voronoi regions only ignores how
        % the objects are grouped in the actual image at time t, since the
        % only Voronoi points added from time t are for new objects
        V_labeled(allmask==0) = 0;
        V_labeled = bwlabel(V_labeled);
        V_merge = zeros(size(V_labeled));
        V_merge_temp = zeros(size(V_labeled));
        
        % V_labeled is tumor image bisected by Voronoi lines into distinct
        % regions. Each of these regions will be marked with the object it
        % belongs to (stored in associated cell)
        associated = cell(1,length(objects));
        used_reg = zeros(1,max(V_labeled(:)));
        for k=1:length(objects)
            pointsX = objects(k).voronoiPointsX;
            pointsY = objects(k).voronoiPointsY;
            associated{k} = zeros(1,length(pointsX));
            tempLL = zeros(size(V_labeled));
            for n = 1:length(pointsX)
                reg = V_labeled(pointsY(n),pointsX(n));
                if(reg ~= 0)
                    used_reg(reg) = used_reg(reg) + 1;
                    associated{k}(n) = reg;
                    tempLL(V_labeled==reg) = 1;
                    if(used_reg(reg)>1)
                        warning('Bad association at region %i: Voronoi point lies on previously assigned region\n',reg);
                    end
                else
                    % this may be happening because the region is too thin
                    fatcut = nonzeros(V_labeled(pointsY(n)+(-1:1),pointsX(n)+(-1:1)));
                    fattercut = nonzeros(V_labeled(pointsY(n)+(-2:2),pointsX(n)+(-2:2)));
                    fattestcut = nonzeros(V_labeled(pointsY(n)+(-3:3),pointsX(n)+(-3:3)));
                    assigned = 0;
                    if(~isempty(fatcut))
                        reg_list = getRegionCandidates(fatcut);
                        reg = reg_list(1);
                        if(used_reg(reg) > 0)
                            % points from object not falling correctly
                            Warning('Ambiguous association at region %i: first choice from fatcut not available\n',reg);
                            % 1) check if an unused region is available
                            for rli = 2:length(reg_list)
                                if(used_reg(reg_list(rli)) == 0)
                                    reg = reg_list(rli);
                                    assigned = 1;
                                    fprintf('Changed association to region %i\n',reg);
                                end
                            end
                        else
                            assigned = 1;
                        end
%                         used_reg(reg) = used_reg(reg) + 1;
%                         associated{k}(n) = reg;
%                         tempLL(V_labeled==reg) = 1;
                        disp(['Minor region ambiguity with ' objs(t).fishID]);
                    end
                    if(~isempty(fattercut) && ~assigned)
                        reg_list = getRegionCandidates(fattercut);
                        reg = reg_list(1);
                        if(used_reg(reg) > 0)
                            % points from object not falling correctly
                            Warning('Ambiguous association at region %i: first choice from fatcut not available\n',reg);
                            % 1) check if an unused region is available
                            for rli = 2:length(reg_list)
                                if(used_reg(reg_list(rli)) == 0)
                                    reg = reg_list(rli);
                                    assigned = 1;
                                    fprintf('Changed association to region %i\n',reg);
                                end
                            end
                        else
                            assigned = 1;
                        end
%                         used_reg(reg) = used_reg(reg) + 1;
%                         associated{k}(n) = reg;
%                         tempLL(V_labeled==reg) = 1;
                        warning(['Medium region ambiguity with ' objs(t).fishID]);
                    end
                    if(~isempty(fattestcut) && ~assigned)
                        reg_list = getRegionCandidates(fattestcut);
                        reg = reg_list(1);
                        if(used_reg(reg) > 0)
                            % points from object not falling correctly
                            Warning('Ambiguous association at region %i: first choice from fatcut not available\n',reg);
                            % 1) check if an unused region is available
                            for rli = 2:length(reg_list)
                                if(used_reg(reg_list(rli)) == 0)
                                    reg = reg_list(rli);
                                    assigned = 1;
                                    fprintf('Changed association to region %i\n',reg);
                                end
                            end
                        else
                            assigned = 1;
                        end
%                         used_reg(reg) = used_reg(reg) + 1;
%                         associated{k}(n) = reg;
%                         tempLL(V_labeled==reg) = 1;
                        warning(['Major region ambiguity with ' objs(t).fishID]);
                    end
                    if(~assigned)
                        % most likely, this object as disappeared from the
                        % image
                        warning(['Apparent disappearing object: ' objs(t).fishID]);
                    else
                        used_reg(reg) = used_reg(reg) + 1;
                        associated{k}(n) = reg;
                        tempLL(V_labeled==reg) = 1;
                    end
                end
            end            
            
            % merge regions
            tempLL = bwmorph(tempLL,'dilate',2);
            tempLL = imfill(tempLL,'holes');
            tempLL = bwmorph(tempLL,'erode',2);
            % avoid overlap
            tempLL(bwmorph(V_merge_temp,'dilate',1)~=0) = 0;
            % add to V_merge labeled with object number
            V_merge_temp = V_merge_temp+k*tempLL;
        end
        % each region should only be assigned once
        if(any(used_reg>1))
            warning(['Bad associations (multiple voronoi points assigned to the same region) may lead to incorrectg segmentation for' objs(end).fishID]);
        end
            
            
        % find unused regions and find associations for them
        unused_reg = find(used_reg==0);
        % 1) if region is attached to a cluster with objects associated
        % with it, find the closest of those objects and assign that
        % 2) if the region cannot be affiliated with any objects, just
        % assign to object with closest Voronoi point
        for u = 1:length(unused_reg)
            reg = unused_reg(u);
            % find if region is attached to an object in original image
            Lreg = return_sub_listL(V_labeled,reg);
            % use Lreg as mask on original image
            Lorig = objs(t).L_labeled;
            Lorig_masked = Lorig;
            Lorig_masked(Lreg==0) = 0;
            % find cluster in Lorig
            clabs = nonzeros(Lorig_masked(:));
            % candidate list of objects (may be pruned later)
            cands = 1:length(objects);
            if(~isempty(clabs))
                % find object to associate with regions
                clab = mode(clabs);
                Lclab = return_sub_listL(Lorig,clab);
                % need to dilate Lclab so that we don't loose good
                % candidates (the issue is images from different times may
                % have a translational shift that introduces error)
                Lclab = bwmorph(Lclab,'dilate',objs(t).threshold);
                % use Lclab as mask to identify candidate objects
                V_cand = V_merge_temp;
                V_cand(Lclab==0) = 0;
                if(any(V_cand(:)>0))
                    cands = nonzeros(unique(V_cand))';
                end
            end
            
            % identify object assignment through closest candidate Voronoi
            rLreg = getCentroids_v1_1(logical(Lreg));
            % sort candidates by distance to region to speed up path
            % distance checks
            cands_dist = zeros(1,length(cands));
            for ik = 1:length(cands)
                k = cands(ik);
                if(~isempty(objects(k).voronoiPointsX))
                    cands_dist(ik) = sqrt(sumsqr(rLreg.Centroid-...
                        [objects(k).voronoiPointsX(end) objects(k).voronoiPointsY(end)]));
                else
                    cands_dist(ik) = inf;
                end
            end
            [~,by_dist] = sort(cands_dist);
            cands = cands(by_dist); %cands are objects, not voronoi points
            kmin = cands(1);
            min_dist = findPathDistance_v1_0([rLreg.Centroid(2) rLreg.Centroid(1)],[objects(kmin).voronoiPointsY(1) objects(kmin).voronoiPointsX(1)],allmask);
            % note the code will repeat the distance check since it hasn't
            % looked at all the voronoi points of the object (don't change
            % the loop without thinking carefully about what is going on!)
            for k = cands
                pointsX = objects(k).voronoiPointsX;
                pointsY = objects(k).voronoiPointsY;
                for n = 1:length(pointsX)
                    dist_n = sqrt(sumsqr(rLreg.Centroid-[pointsX(n) pointsY(n)]));
                    if(dist_n <= min_dist)
                        % check path distance
                        pdist_n = findPathDistance_v1_0([rLreg.Centroid(2) rLreg.Centroid(1)],[pointsY(n) pointsX(n)],allmask);
                        %                         [pdist_n,anchors,~] = findPathDistance_v1_0([rLreg.Centroid(2) rLreg.Centroid(1)],[pointsY(n) pointsX(n)],allmask);
                        %                         if(any(weights<0))
                        %                             error('Negative Weights');
                        %                         end
                        if(pdist_n < min_dist)
                            %                             imshow(logical(V_labeled));
                            %                             hold on;
                            %                             plot([rLreg.Centroid(1); anchors(:,2); objects(k).voronoiPointsX(n)],[rLreg.Centroid(2); anchors(:,1); objects(k).voronoiPointsY(n)],'+-',objects(k).voronoiPointsX(:),objects(k).voronoiPointsY(:),'*');
                            %                             hold off;
                            %                             pause(0.1);
                            min_dist = pdist_n;
                            kmin = k;
                        end
                    end
                end
            end
            
            % add to associated cell
            associated{kmin} = [associated{kmin} reg];
        end
        
        % create V_merge
        for k=1:length(objects)
            
            % make empty matrix with 0's
            tempLL = zeros(size(V_labeled));
            
            % fill in 1's in the regions which need to be joined
            for g=1:length(associated{k})
                if(associated{k}(g)~=0)
                    tempLL(V_labeled==associated{k}(g))=1;
                else
                    % this bug happens if a voronoi point is in the empty
                    % region for some reason. Could be assosiated with a
                    % disappearing tumor warning
                    warning('Voronoi point in empty space at %s',objs(t).fishID);
                end
            end
            
            % merge regions
            tempLL = bwmorph(tempLL,'dilate',2);
            tempLL = imfill(tempLL,'holes');
            tempLL = bwmorph(tempLL,'erode',2);
            % avoid overlap
            tempLL(bwmorph(V_merge,'dilate',1)~=0) = 0;
            % add to V_merge labeled with object number
            V_merge = V_merge+k*tempLL;
        end
        
        %
        % %             figure(2)
        % %             imshow(objs(t).L_t1*63/nmets(t),colormap('hot'));
        % %             Bt = bwboundaries(tempLL);
        % %             hold on
        % %             plot(pointsX,pointsY,'*')
        % %             for kk=1:length(Bt)
        % %                 b=Bt{kk};
        % %                 plot(b(:,2),b(:,1),'-g')
        % %             end
        % %             pause(0.1);
        
        L_divideds{t} = V_merge;
        L_divideds{t}(objs(t).L_t1==0)=0;
        
    end

    
    rC = getCentroids_v1_1(L_divideds{t-1}); % get CM's of objects in previous time point
    % find out which clusters are actually here
    rClist = nonzeros(unique(L_divideds{t-1}))';
    % mark L_divs(t) with dots from tumors at L_divs(t-1) if no other tumor
    % is present within threshold of centroid
    L_save = L_divideds{t};
    L_list = nonzeros(unique(L_save));
    if(isempty(L_list))
        L_max = 0;
    else
        L_max = max(L_list);
    end
    for k = rClist
        % This code assumes the indices of rClist correspond to the object
        % number the clusters belong to - this works because V_merge is
        % labeled with the object numbers
        
        L_temp = L_save;
        L_temp(V_merge~=k) = 0;
        if(~any(L_temp(:)))
            % cluster is not persistent; add a pointer
            rcent = rC(k).Centroid;
            % use voronoi point closest to rcent
            if(exist('objects','var') && ~isempty(objects(k).voronoiPointsX))
                vpoints = [objects(k).voronoiPointsX objects(k).voronoiPointsY];
                dist_sqr_vp = sum((vpoints - rcent).^2,2);
                [~,vpi] = min(dist_sqr_vp);
                cent = vpoints(vpi,:);
            else
                cent = round(rcent);
            end
            L_divideds{t}(cent(2),cent(1)) = L_max + 1;
            L_max = L_max + 1;
        end
    end
 
    
    % fill cross_index matrix
    % need to make sure clusters grouped before same as clusters below...
    rL = getCentroids_v1_1(L_divideds{t});   % get CM's of blobs in this time point
    % find out which clusters are actually here
    rLlist = nonzeros(unique(L_divideds{t}))';
    
    if(length(rLlist)<length(rClist))
        % bad job finding metastasis that have disappeared
        error(['Failed cluster allocation with ' objs(end).fishID]);
    end
    
    cmdists = zeros(length(rClist),length(rLlist));
    
    % get distances between clusters
    for k = 1:length(rClist)
        for n = 1:length(rLlist)
            cent1 = rC(rClist(k)).Centroid;
            cent2 = rL(rLlist(n)).Centroid;
            cmdists(k,n) = sqrt(sumsqr(cent1-cent2));
        end
    end
    priorities = zeros(size(cmdists));
    for k = 1:length(rClist)
        [~,priorities(k,:)] = sort(cmdists(k,:));
    end
    % get arrangement with no conflicts that (somewhat) minimizes
    % distance
    if(~isempty(priorities))
        cross_links{t-1} = findConfig_v1_3(priorities,ones(1,length(rClist)),cmdists);
    end
    % the unassigned clusters (new metastasis)
    unlinked = setdiff(1:length(rLlist),cross_links{t-1});
    % rename clusters in L_divided for consistency
    newname = zeros(1,max(rLlist));

    for k = 1:length(rClist)
        newname(rLlist(cross_links{t-1}(k))) = rClist(k);
    end
    if(isempty(rClist))
        rCmax = 0;
    else
        rCmax = max(rClist);
    end
    
    for nn = 1:length(unlinked)
        newname(rLlist(unlinked(nn))) = rCmax + nn;
    end
    
    % make cross_index; since the L_divideds are being renamed, the
    % cross_index is actually trivial (same number in all rows). However, I
    % will still calculate it for 1) robustness of code in case I change
    % something in the future and 2) as a check to see renaming went well
    for tt = 1:t-1
        for ii = 1:length(index{tt})
            ci1 = cross_index{tt}(ii,t-tt);
            k = find(rClist==ci1,1);
            cross_index{tt}(ii,t-tt+1) =  newname(rLlist(cross_links{t-1}(k)));
        end
    end
    % need to account for new clusters
    %     disp(nmets);
    %     disp(rLlist(cross_links{t-1}));
    %     disp(rLlist(unlinked));
    %     disp(newname);
%     imshow(L_divideds{t}*36/nmets(t),colormap('hot'));
%     pause(0.1);
    if(nmets(t)~=(length(unlinked)+length(cross_links{t-1})))
        error(['Segmentation Error 1 with ' objs(end).fishID]);
    end
    cross_index{t}(:,1) = newname(rLlist(unlinked));

    % actually renaming:
    znewname = [0 newname];
    L_divideds{t} = znewname(L_divideds{t}+1);
    
    rLmax2 = length(regionprops(L_divideds{t}));
    if(max(newname)~=rLmax2)
        error(['Segmentation Error 2 with ' objs(end).fishID]);
    end
    
    %     if(nmets(t)>1)
    %         figure(1);
    %         rows = 150:450;
    %         cols = 550:1000;
    %         imshow(double(L_divideds{t}(rows,cols))*63/nmets(t),colormap('hot'));
    % %         hold on;
    % %         bounds = bwboundaries(objs(t).L_t1);
    % %         for c = 1:length(bounds)
    % %             b = bounds{c};
    % %             plot(b(:,2)-cols(1)+1,b(:,1)-rows(1)+1,'w-','LineWidth',1.25);
    % %         end
    %         pause(0.1);
    %     end
    
    % if V_merge is empty (meaning L_D7 did not need to get divided) then funtion should just return one big region for
    % entire image
    if isempty(V_merge)==1
        V_merges{t-1} = ones(size(objs(t).L_t1));
    else
        V_merges{t-1} = V_merge;
    end
    
    V_merge=[]; % V_merge empty when voronoi not used
    
end


end

function regs = getRegionCandidates(fatcut)

cands = unique(fatcut);
pop = zeros(size(cands));

for i = 1:length(cands)
    c = cands(i);
    pop(i) = sum(fatcut==c);
end

[~,ord] = sort(pop);
regs = cands(flipud(ord));
end



