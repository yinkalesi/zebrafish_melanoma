function [rev_tum,rev_Ldivs] = getL_Revised_v1_4(tumors,images,old_rev)
%% getL_Revised
%  version 1.4
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 2/20/20

%% Version History
%  1.0: this lets a user manually screen tumors 
%  1.1: tumor struct also revised
%  1.2: better interface!
%  1.3: can take existing rev_tum as input so you have the option to use
%  previous revisions to speed up new revisions (check to make sure things
%  are compatible)
%  1.4: making composite image instead of using imshowpair
%  3/11/21: changed images(i,n).bf to images(i,n).bf(:,:,1) to accomodate
%  Silja's data

% 1) Iterate through each Ldiv
% 2) For each region in final Ldiv, show the tumor history and verify if it
% is a real tumor

if(exist('old_rev','var'))
    record_present = true;
else
    record_present = false;
    old_rev = struct();
end
rev_Ldivs = cell(size(images));

% initialize rev_tum
rev_tum = tumors;
% tumors = struct('sizes',cell(1,nfsh),'locs',cell(1,nfsh),...
%     'reorg',cell(1,nfsh),'count',cell(1,nfsh),...
%     'distances',cell(1,nfsh),'length',cell(1,nfsh),...
%     'width',cell(1,nfsh),'perimeter',cell(1,nfsh),...
%     'intensity',cell(1,nfsh),'intensity_std',cell(1,nfsh),...
%     'convex_solidity',cell(1,nfsh),'gender',cell(1,nfsh),...
%     'absorbance',cell(1,nfsh),'body_brightness',cell(1,nfsh),...
%     'absorbance2',cell(1,nfsh),'body_brightness2',cell(1,nfsh));

screen_size = get(0,'ScreenSize');
fig1 = figure;
fig1.Position = screen_size;
fig2 = figure;
fig2.Position = screen_size;
for n = 1:size(images,2)
    % find last time these is an image
    if(~isempty(images(1,n).count))
        tf = images(1,n).count;
    else
        tf = 0;
    end
    
    % iterate through all regions
    sizes = tumors(n).sizes;
    locs = tumors(n).locs;
    nr = size(locs,1);
    remove_list = zeros(1,nr);
    rmc = 0;
    figure(fig1);
    picx = min(size(images(1,n).Ldiv,2),400);
    picy = min(size(images(1,n).Ldiv,1),300);
    for r = 1:nr
        clf(fig1,'reset');
        % iterate through all clusters and store verified clusters in Lrev
        t0 = 1;
        t1 = tf;
        npic_max = 12;
        if(tf>npic_max)
            % consider showing subselection of pics
            t0 = find(sizes(r,:)>0,1,'first');
            if(isempty(t0))
                t0 = 1;
            elseif(t0>1)
                t0 = t0-1; % first frame shown will now have metastasis region
            end
            t1 = min(tf,t0-1+npic_max);
        end
        nti = t1-t0+1;
        nrow = ceil((nti-0.5)/3);
        ncol = ceil((nti-0.5)/nrow);
        for i = t0:t1
            reg = return_sub_listL(images(i,n).Ldiv,r);
            if(~any(reg(:)))
                subplot(nrow,ncol,i-t0+1);
                picbox = [max(1,round(locs(r,2*tf-1:2*tf)-[picx picy]*0.5)) picx picy];
                scale = 128/double(max(max(imcrop(images(i,n).fimg,picbox))));
                comp = uint8(zeros(size(images(i,n).rgb)));
                comp(:,:,1) = scale*images(i,n).auto+0.5*images(i,n).bf(:,:,1);
                comp(:,:,2) = scale*images(i,n).fimg+0.5*images(i,n).bf(:,:,1);
                comp(:,:,3) = 0.5*images(i,n).bf(:,:,1);
                imshow(imcrop(comp,picbox));
                %                 imshowpair(imcrop(images(i,n).fimg,picbox),imcrop(images(i,n).auto,picbox),'ColorChannels','red-cyan');
                %                	imshow(false([picy picx]));
                title(sprintf('f %i, t %i, r %i not present',n,i,r));
            else
                reg = bwmorph(reg,'dilate',2);
                max_int = max(images(i,n).fimg(logical(reg(:))));
                breg = bwboundaries(reg);
                subplot(nrow,ncol,i-t0+1);
                picbox = [max(1,round(locs(r,2*i-1:2*i)-[picx picy]*0.5)) picx picy];
                scale = 128/double(max(max(imcrop(images(i,n).fimg,picbox))));
                comp = uint8(zeros(size(images(i,n).rgb)));
                comp(:,:,1) = scale*images(i,n).auto+0.5*images(i,n).bf(:,:,1);
                comp(:,:,2) = scale*images(i,n).fimg+0.5*images(i,n).bf(:,:,1);
                comp(:,:,3) = 0.5*images(i,n).bf(:,:,1);
                imshow(imcrop(comp,picbox));
                %                 imshowpair(imcrop(images(i,n).fimg,picbox),imcrop(images(i,n).auto,picbox),'ColorChannels','red-cyan');
                hold on;
                for c=1:length(breg)
                    b=breg{c};
                    plot(b(:,2)-picbox(1)+1,b(:,1)-picbox(2)+1,'y-','LineWidth',0.5)
                end
                text(locs(r,2*i-1)-picbox(1)+1,locs(r,2*i)-picbox(2)+1,['\leftarrow' num2str(r) ': ' num2str(sizes(r,i))],'Color',[1 1 1]);
                hold off;
                title(sprintf('f %i, t %i, r %i of %i, (area,I_{max})=(%i,%i)',n,i,r,nr,sizes(r,i),max_int));
            end
        end
        fig1.Position = screen_size;
        if(~record_present)
            % pop up dialog box
            qres = questdlg('Is this a real tumor?','Signal Verification','Yes','No','Reject Remaining','Yes');
            switch qres
                case 'Yes'
                    % add region to Lrev
                case 'No'
                    % don't add region
                    rmc = rmc + 1;
                    remove_list(rmc) = r;
                case 'Reject Remaining'
                    % remove all remaining regions
                    remove_list(rmc+1:rmc+1+nr-r) = r:nr;
                    rmc = rmc+1+nr-r;
                    break;
            end
        else
            % pop up dialog box
            qres = questdlg('Is this a real tumor?','Signal Verification','Yes','No','Use Record','Yes');
            switch qres
                case 'Yes'
                    % add region to Lrev
                case 'No'
                    % don't add region
                    rmc = rmc + 1;
                    remove_list(rmc) = r;
                case 'Use Record'
                    % add clusters in record with id >= r to remove list
                    old_rem_list = old_rev(n).remove_list;
                    start_r = find(old_rem_list>=r,1,'first');
                    if(~isempty(start_r))
                        n_add = length(old_rem_list)-start_r+1;
                        remove_list(rmc+1:rmc+n_add) = old_rem_list(start_r:end);
                        rmc = rmc+n_add;
                    end
                    break;
            end
        end
    end
    remove_list = remove_list(1:rmc);
    rev_tum(n).remove_list = remove_list;
%     close;
    
    % get not_remove list
    not_removed = zeros(1,nr-rmc);
    nri = 0;
    % first chunk
    if(rmc>0)
        not_removed(1:remove_list(1)-1) = 1:remove_list(1)-1;
        nri = remove_list(1)-1;
    end
    % intermediate chunks
    for rli = 2:rmc
        rldiff = remove_list(rli)-remove_list(rli-1);
        not_removed(nri+1:nri+rldiff-1) = remove_list(rli-1)+1:remove_list(rli)-1;
        nri = nri+rldiff-1;
    end
    % final chunk
    if(rmc>0)
        not_removed(nri+1:end) = remove_list(rmc)+1:nr;
    else
        not_removed(nri+1:end) = 1:nr;
    end
        
    % make revised Ldivs and tumor struct at each time point
    for i = 1:tf
        Ldiv = images(i,n).Ldiv;
        Lrev = zeros(size(Ldiv),'uint8');
        for nri = 1:length(not_removed)
            Lrev(Ldiv==not_removed(nri)) = nri;
        end
%         for p = 1:length(Lrev(:))
%             ireg = find(not_removed==Ldiv(p));
%             if(~isempty(ireg))
%                 % not to be removed
%                 Lrev(p) = ireg;
%             end
%         end
        rev_Ldivs{i,n} = Lrev;
    end
    
    % revise appropriate field in tumors based on remove_list
    fieldsToChange = {'sizes','locs','length','width','perimeter',...
        'intensity','intensity_std','convex_solidity','absorbance',...
        'absorbance2'};
    
    for ftci = 1:length(fieldsToChange)
        field = fieldsToChange{ftci};
        data = rev_tum(n).(field);
        % do update
        new_data = data(not_removed,:);
        rev_tum(n).(field) = new_data;
    end
    
    % need to change tumor count by number removed
    for i = 1:tf
        for r = remove_list
            Ldiv = images(i,n).Ldiv;
            if(any(Ldiv(:)==r))
                rev_tum(n).count(i) = rev_tum(n).count(i)-1;
            end
        end
    end
    
    figure(fig2);
    clf(fig2,'reset');
    % plot revised image
    sum_img = zeros(size(images(1,n).Ldiv));
    for i = 1:tf
        sum_img = sum_img+logical(rev_Ldivs{i,n});
    end
    if(tf > 0)
        scale = 0.75;
        imshow(imresize(sum_img*(63/tf),scale),hot);
        hold on;
        for r = 1:length(not_removed)
            % find position of identified tumor
            mloc = rev_tum(n).locs(r,2*tf-1:2*tf);
            % print fid number
            text(mloc(1)*scale,mloc(2)*scale,['\leftarrow' num2str(r)],'Color',[1 1 1]);
        end
        hold off;
        title(sprintf('f %i, nr %i',n,length(not_removed)));
        pause(1);
    end
end
