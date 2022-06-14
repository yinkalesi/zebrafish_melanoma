%% getScreenedTumorDistributions_v4_0
%  Version 4.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 8/21/20
% take data from several summary files and create one distribution
%% Version History
%  Adds a screening step to combineTumorSummarys (sic) function
%  1.1: use average intensity data as well
%  1.2: added option the remove certain fish and a few minor mods
%  1.3: makes sure can combine tumor files...
%  1.4: separate by gender - pick only male or female fish
%  1.5: handle split up mat files
%  1.6: changing names to handle 'revised' data, not using totals
%  1.7: using absorbance2
%  2.0: using size_range to select only fish with day1 tumor area of a
%  certain size; set is_disappeared to always be false (user revision
%  process should obviate this need)
%  3.0: return which are male and female
%  3.1: adding a new count+averaging method whis is nfish+cured fish
%  producing extra file to identify primaries over time
%  3.2: adding check to detect and remove artifacts (objects that appear on
%  day 1 but aren't there later)
%  3.3: raised cured fish criteria to x<=1000
%  4.0: change so size_range works on an individual tumor basis rather than
%  whole fish
%  4.1: raised cured fish criteria to x<=1000 (updated version of 3.3);
%  fixed calculating active count so everything is a whole number + we
%  account for when batches of fish were measured

function [dist,x,dist2,x2,fish_desc] = getScreenedTumorDistributions_v4_1(tumor_files,times,name,removal_list,gender_id,size_range,min_survival_time)

summary_files(length(tumor_files)) = struct('tumors',struct());
ntums = zeros(1,length(tumor_files));
nprims = zeros(1,length(tumor_files));
nfish = zeros(1,length(tumor_files));
t = [];
for itf = 1:length(tumor_files)
    temp = load([tumor_files{itf} '_tumors_revised.mat']);
    summary_files(itf).tumors = temp.rev_tum;
    
    for itt = 1:length(summary_files(itf).tumors)
        if(any(summary_files(itf).tumors(itt).count))
            ntums(itf) = ntums(itf)+summary_files(itf).tumors(itt).count(end);
            nprims(itf) = nprims(itf)+summary_files(itf).tumors(itt).count(1);
            nfish(itf) = nfish(itf)+logical(summary_files(itf).tumors(itt).count(1));
        else
            fprintf('File %i, fish %i has no data\n',itf,itt);
        end
    end
    t = unique([t times{itf}]);
    
%     temp = load([tumor_files{itf} '_totals.mat']);
%     summary_files(itf).totals = temp.totals;
%     ntums(itf) = size(summary_files(itf).totals.distrib,1);
%     nprims(itf) = sum(summary_files(itf).totals.primary_counts);
%     nfish(itf) = sum(summary_files(itf).totals.fish_counts);
%     t = unique([t times{itf}]);
end

max_time_ind = length(t);
% data sets may not use the same time increments - so need to find the
% right time point indices for each data set
time_indices = cell(1,length(tumor_files));
meas_at_t = false(length(tumor_files),max_time_ind);
for itf = 1:length(tumor_files)
    tref = zeros(size(times{itf}));
    itref = 1;
    for it = 1:max_time_ind
        if(t(it)==times{itf}(itref))
            tref(itref) = it;
            itref = itref+1;
            meas_at_t(itf,it) = true;
            if(itref > length(tref))
                break;
            end
        end
    end
    time_indices{itf} = tref;
end

ntums_total = sum(ntums);
dist = zeros(ntums_total,max_time_ind);
x = zeros(ntums_total,1);
yz = zeros(ntums_total,1);
nprims_total = sum(nprims);
dist2 = zeros(nprims_total,max_time_ind);
x2 = zeros(nprims_total,1);
yz2 = zeros(nprims_total,1);
nfish_total = sum(nfish);
nlargest_sizes = nfish_total*max_time_ind;
x3 = zeros(nlargest_sizes,1);
yz3 = zeros(nlargest_sizes,1);
dist3 = zeros(nlargest_sizes,max_time_ind);
fish_counts = zeros(1,max_time_ind);
cured_fish = zeros(length(tumor_files),max_time_ind);
cured_prim = zeros(length(tumor_files),max_time_ind);
prim_counts = zeros(1,max_time_ind);
tumor_counts = zeros(1,max_time_ind);
itum = 0;
iprim = 0;
ilarge = 0;
% want to store size and intensity values in one place
sizes = zeros(ntums_total,max_time_ind);
attenuations = zeros(ntums_total,max_time_ind);
intensities = zeros(ntums_total,max_time_ind);
store_labels = zeros(ntums_total,3); %fid tfid isLargest?
% fish_disc: tumor_file,fid,gender,t1_area,tlast
fish_desc = zeros(nfish_total,5);
nfish_sel = 0;
nprim_sel = 0;
i_store = 0;
nfish_used = zeros(length(tumor_files),1);
for itf = 1:length(tumor_files)
    summary = summary_files(itf);
    toRemove = removal_list{itf};
    tref = time_indices{itf};
    % iterate through tumors and screen individually
    % only analyze one gender...
    nfshi = length(summary.tumors);
    % get file descriptor
    bslashpos = find(tumor_files{itf}=='\',1,'last');
    if(isempty(bslashpos))
        bslashpos = 0;
    end
    tm_desc = tumor_files{itf}(bslashpos+1:end);
    for fid = 1:nfshi
        tm_sizef = summary.tumors(fid).sizes;
        tm_intf = summary.tumors(fid).intensity;
        
        if(isfield(summary.tumors(fid),'absorbance2'))
            tm_absf = summary.tumors(fid).absorbance2;
        else
            tm_absf = zeros(size(tm_intf));
        end
        % screen based on if considered an 'artifact'
        % tumors that appear on time 1 but not time 2 need to be screened
        % out
        all_tum_disappeared = tm_sizef(:,:)==1;
        all_present_at_day1 = tm_sizef(:,1)~=0;
        if(size(all_tum_disappeared,2)>1)
            all_not_present_enough = all_tum_disappeared(:,2)==1;
        else
            all_not_present_enough = false(size(all_tum_disappeared(:,1)));
        end
        
        possible_artifact = all_not_present_enough&all_present_at_day1;
        all_possible_artifacts = all(possible_artifact);
        
        if(any(possible_artifact))
            fprintf('Disappeared signals at file %s, fid %i, pos %s were removed\n',...
                tm_desc,fid,num2str(find(possible_artifact)'));
        end
        
        tm_sizef = tm_sizef(~possible_artifact,:);
        tm_intf = tm_intf(~possible_artifact,:);
        tm_absf = tm_absf(~possible_artifact,:);
        
        % screen based on if tumors exist within size range
        tum_in_size_range = (tm_sizef(:,1)>=size_range(1))&(tm_sizef(:,1)<=size_range(2));
        if(any(~tum_in_size_range))
            fprintf('Out of range signals at file %s, fid %i, adj pos %s were removed\n',...
                tm_desc,fid,num2str(find(~tum_in_size_range)'));
        end
        tm_sizef = tm_sizef(tum_in_size_range,:);
        tm_intf = tm_intf(tum_in_size_range,:);
        tm_absf = tm_absf(tum_in_size_range,:);
        
        [ntumsf,tlast] = size(tm_sizef);
        tm_count = sum(tm_sizef>0,1);
        
        f_gender = summary.tumors(fid).gender;
        if(~isempty(gender_id) && ~isempty(f_gender))
            is_right_gender = any(f_gender==gender_id);
        else
            is_right_gender = true;
        end
        % ignore fish that died the first day
        % check if tumor area on first day is within size range
        if(tlast)
            t1_area = sum(tm_sizef(:,1));
        else
            t1_area = 0;
        end
        if(any(tum_in_size_range) && ~isempty(tm_sizef))
            is_in_range = true;
        else
            is_in_range = false;
        end
        survived_long_enough = tlast>=min_survival_time;
        
        if(survived_long_enough && ~any(toRemove==fid) && is_right_gender && is_in_range && ~all_possible_artifacts)
            
            fish_counts(tref(1:tlast)) = fish_counts(tref(1:tlast)) + 1;
            prim_counts(tref(1:tlast)) = prim_counts(tref(1:tlast)) + tm_count(1);
            % check for cured fish
            if(all(tm_sizef(:,tlast)<=1000))
                cured_fish(itf,tref(tlast)) = cured_fish(itf,tref(tlast)) + 1;
                cured_prim(itf,tref(tlast)) = cured_prim(itf,tref(tlast)) + tm_count(1);
                %                 nprims_active(tref(tlast+1:end)) = nprims_active(tref(tlast+1:end)) + tm_count(1);
            end
            nfish_sel = nfish_sel + 1;
            nprim_sel = nprim_sel + tm_count(1);
            fish_desc(nfish_sel,1:5) = [itf,fid,f_gender,t1_area,tlast];
            nfish_used(itf) = nfish_used(itf)+1;
            % find largest tumor (will be designated a primary)
            [~,tfi_largest] = max(tm_sizef(:,1));
            for it = 1:tlast
                ilarge = ilarge+1;
                x3(ilarge) = tm_sizef(tfi_largest,it);
                % calculate attenuation due to pigmentation
                A_pig = tm_absf(tfi_largest,it)-tm_absf(tfi_largest,1);
                k_atten = getAttenuation(A_pig);
                yz3(ilarge) = x3(ilarge)*tm_intf(tfi_largest,it)/k_atten;
                dist3(ilarge,tref(it)) = 1;
                store_labels(i_store+tfi_largest,3) = 1;
            end
            for tfid = 1:ntumsf
                % find time tumor appeared
                itstart = find(tm_sizef(tfid,:)==0,1,'last')+1;
                if(isempty(itstart))
                    itstart = 1;
                end
                for it = itstart:tlast
                    itum = itum+1;
                    x(itum) = tm_sizef(tfid,it);
                    % calculate attenuation due to pigmentation
                    A_pig = tm_absf(tfid,it)-tm_absf(tfid,itstart);
                    k_atten = getAttenuation(A_pig);
                    yz(itum) = x(itum)*tm_intf(tfid,it)/k_atten;
                    sizes(i_store+tfid,tref(it)) = x(itum);
                    attenuations(i_store+tfid,tref(it)) = k_atten;
                    intensities(i_store+tfid,tref(it)) = yz(itum);
                    store_labels(i_store+tfid,1:2) = [fid tfid];
                    dist(itum,tref(it)) = 1;
                    tumor_counts(tref(it)) = tumor_counts(tref(it)) + 1;
                end
                % add to primary list if this existed on Day 1
                if(itstart==1)
                    for it = 1:tlast
                        iprim = iprim+1;
                        x2(iprim) = tm_sizef(tfid,it);
                        % calculate attenuation due to pigmentation
                        A_pig = tm_absf(tfid,it)-tm_absf(tfid,1);
                        k_atten = getAttenuation(A_pig);
                        yz2(iprim) = x2(iprim)*tm_intf(tfid,it)/k_atten;
                        dist2(iprim,tref(it)) = 1;
                        %                             primary_counts(tref(it)) = primary_counts(tref(it)) + 1;
                    end
                end
            end
            i_store = i_store+ntumsf;
        else
            if(~survived_long_enough)
                fprintf('File %s, Fish %i did not survive at least %i days and was ignored\n',tm_desc,fid,min_survival_time);
            elseif(~is_right_gender)
                fprintf('File %s, Fish %i was ignored due to gender restriction\n',tm_desc,fid);
            elseif(all_possible_artifacts)
                fprintf('File %s, Fish %i was ignored due to inconsistent signals\n',tm_desc,fid);
            elseif(~is_in_range)
                fprintf('File %s, Fish %i was ignored due to size restrictions\n',tm_desc,fid);
            else
                fprintf('File %s, Fish %i was ignored due to removal list\n',tm_desc,fid);
            end
        end
    end
end

% calculate active count based on adding back cured fish

cum_cured_fish_batch = [zeros(size(cured_fish,1),1) cumsum(cured_fish(:,1:end-1),2)];
cum_cured_prim_batch = [zeros(size(cured_prim,1),1) cumsum(cured_prim(:,1:end-1),2)];

% calculate eff_meas which is meas_at_t except times we would have measured
% the fish if they had not been cured are also marked as true. We measure
% the fish every odd day so every odd day after the last measurement of
% each batch is marked as true
odd_times = mod(t,2)==1;
eff_meas = meas_at_t;
for itf = 1:size(meas_at_t,1)
    last_meas_i = find(meas_at_t(itf,:),1,'last');
    eff_meas(itf,odd_times&((1:length(t))>last_meas_i)) = true;
end

total_cured_fish = sum(cum_cured_fish_batch.*eff_meas,1);
total_cured_prim = sum(cum_cured_prim_batch.*eff_meas,1);

% total_cured_fish = sum(cured_fish,1);
% total_cured_prims = sum(cured_prim,1);
% cum_cured_fish = [0 cumsum(total_cured_fish(1:end-1),2)];
% cum_cured_prim = [0 cumsum(total_cured_prims(1:end-1),2)];
% % cum_cured_fish = [zeros(size(cured_fish,1),1) cumsum(cured_fish(:,1:end-1),2)];
% % cum_cured_prim = [zeros(size(cured_prim,1),1) cumsum(cured_prim(:,1:end-1),2)];
% active_fish = nfish_sel*fish_counts./(nfish_sel-cum_cured_fish);
% active_prim = nprim_sel*prim_counts./(nprim_sel-cum_cured_prim);
% % for itf = 1:length(tumor_files)
% %     tref = time_indices{itf};
% %     active_fish(tref) = active_fish(tref)+cum_cured_fish(itf,tref);
% % end

active_fish = fish_counts+total_cured_fish;
active_prim = prim_counts+total_cured_prim;

disp(['   Fish Counts: ' num2str(fish_counts)]);
disp(['Primary Counts: ' num2str(prim_counts)]);
disp(['  Tumor Counts: ' num2str(tumor_counts)]);
disp(['    Cured Fish: ' num2str(total_cured_fish)]);
disp([' Cured Primary: ' num2str(total_cured_prim)]);
disp(['   Active Fish: ' num2str(active_fish)]);
disp(['Active Primary: ' num2str(active_prim)]);

% truncate to correct size
sizes = sizes(1:i_store,:);
attenuations = attenuations(1:i_store,:);
intensities = intensities(1:i_store,:);
store_labels = store_labels(1:i_store,:); % [fid tfid is_prim]
fish_desc = fish_desc(1:nfish_sel,:); % [itf,fid,f_gender,t1_area,tlast];
% % select largest tumors from stored values
% select_largest = find(store_labels(:,3)==1);
% sizes_large = sizes(select_largest,:);
% atten_large = attenuations(select_largest,:);
% inten_large = intensities(select_largest,:);
% % plot some results
% figure;
% plot(t,sizes_large);
% figure;
% subplot(1,2,1);
% plot(t,inten_large);
% subplot(1,2,2);
% plot(t,inten_large.*atten_large);
% figure;
% plot(t,atten_large);

dotloc = find(name=='.',1); % need in order to make 'name' a prefix
% store counts
dlmwrite([name(1:dotloc-1) '_counts.txt'],[fish_counts; prim_counts;...
    tumor_counts; total_cured_fish; active_fish; active_prim]);
% store size and intensity trajectories
dlmwrite([name(1:dotloc-1) '_sizes.txt'],sizes);
dlmwrite([name(1:dotloc-1) '_attenuation_values.txt'],attenuations);
dlmwrite([name(1:dotloc-1) '_intensities2.txt'],intensities);
dlmwrite([name(1:dotloc-1) '_labels.txt'],store_labels);
dlmwrite([name(1:dotloc-1) '_desc.txt'],fish_desc);

[~,ord] = sort(x(1:itum));
[~,ord2] = sort(x2(1:iprim));
[~,ord3] = sort(x3(1:ilarge));
[~,ord4] = sort(yz(1:itum));
[~,ord5] = sort(yz2(1:iprim));
[~,ord6] = sort(yz3(1:ilarge));
temp_x = x(ord);
temp_x2 = x2(ord2);
temp_x3 = x3(ord3);
temp_yz = yz(ord4);
temp_yz2 = yz2(ord5);
temp_yz3 = yz3(ord6);
temp_dist = dist(ord,:);
temp_dist2 = dist2(ord2,:);
temp_dist3 = dist3(ord3,:);
temp_dist4 = dist(ord4,:);
temp_dist5 = dist2(ord5,:);
temp_dist6 = dist3(ord6,:);

% need to combine numbers if they are the same
prev = -1;
pos = 0;
x = zeros(size(temp_x));
dist = zeros(size(temp_dist));
for i = 1:length(temp_x)
    if(temp_x(i)==prev)
        dist(pos,:) = dist(pos,:)+temp_dist(i,:);
    else
        prev = temp_x(i);
        pos = pos+1;
        x(pos) = temp_x(i);
        dist(pos,:) = dist(pos,:)+temp_dist(i,:);
    end
end
x = x(1:pos);
% turn 1 into zero
if(pos>0 && x(1)==1)
    x(1)=0;
end
dist = dist(1:pos,:);

prev = -1;
pos = 0;
yz = zeros(size(temp_yz));
dist4 = zeros(size(temp_dist4));
for i = 1:length(temp_yz)
    if(temp_yz(i)==prev)
        dist4(pos,:) = dist4(pos,:)+temp_dist4(i,:);
    else
        prev = temp_yz(i);
        pos = pos+1;
        yz(pos) = temp_yz(i);
        dist4(pos,:) = dist4(pos,:)+temp_dist4(i,:);
    end
end
yz = yz(1:pos);
dist4 = dist4(1:pos,:);

saveDistribution([name(1:dotloc-1) '_area_vs_count' name(dotloc:end)],dist,x,t);
saveDistribution([name(1:dotloc-1) '_intensity2_vs_count' name(dotloc:end)],dist4,yz,t);

% make normalized distributions
% 1) normalize by number of fish at time t
dist_nbf = getNorm(dist,fish_counts);
dist4_nbf = getNorm(dist4,fish_counts);
saveDistribution([name(1:dotloc-1) '_area_vs_count_over_fish' name(dotloc:end)],dist_nbf,x,t);
saveDistribution([name(1:dotloc-1) '_intensity2_vs_count_over_fish' name(dotloc:end)],dist4_nbf,yz,t);
% 2) normalize by number of Day 1 tumors in fish at time t
dist_nbp = getNorm(dist,prim_counts);
dist4_nbp = getNorm(dist4,prim_counts);
saveDistribution([name(1:dotloc-1) '_area_vs_count_over_prim' name(dotloc:end)],dist_nbp,x,t);
saveDistribution([name(1:dotloc-1) '_intensity2_vs_count_over_prim' name(dotloc:end)],dist4_nbp,yz,t);
% 3) normalize by number of tumors at time t
dist_nbt = getNorm(dist,tumor_counts);
dist4_nbt = getNorm(dist4,tumor_counts);
saveDistribution([name(1:dotloc-1) '_area_vs_count_over_tum' name(dotloc:end)],dist_nbt,x,t);
saveDistribution([name(1:dotloc-1) '_intensity2_vs_count_over_tum' name(dotloc:end)],dist4_nbt,yz,t);
% 4) normalize by number of 'active fish': (nfish+number that were cured in
% the past)
dist_nba = getNorm(dist,active_fish);
dist4_nba = getNorm(dist4,active_fish);
saveDistribution([name(1:dotloc-1) '_area_vs_count_over_active' name(dotloc:end)],dist_nba,x,t);
saveDistribution([name(1:dotloc-1) '_intensity2_vs_count_over_active' name(dotloc:end)],dist4_nba,yz,t);
% 5) normalize by number of 'active primaries': (nprim+nprim in fish that 
% were cured in the past)
dist_nbap = getNorm(dist,active_prim);
dist4_nbap = getNorm(dist4,active_prim);
saveDistribution([name(1:dotloc-1) '_area_vs_count_over_actprim' name(dotloc:end)],dist_nbap,x,t);
saveDistribution([name(1:dotloc-1) '_intensity2_vs_count_over_actprim' name(dotloc:end)],dist4_nbap,yz,t);

% need to get primary dist
prev = -1;
pos = 0;
x2 = zeros(size(temp_x2));
dist2 = zeros(size(temp_dist2));
for i = 1:length(temp_x2)
    if(temp_x2(i)==prev)
        dist2(pos,:) = dist2(pos,:)+temp_dist2(i,:);
    else
        prev = temp_x2(i);
        pos = pos+1;
        x2(pos) = temp_x2(i);
        dist2(pos,:) = dist2(pos,:)+temp_dist2(i,:);
    end
end
x2 = x2(1:pos);
% turn 1 into zero
if(x2(1)==1)
    x2(1)=0;
end
dist2 = dist2(1:pos,:);

prev = -1;
pos = 0;
yz2 = zeros(size(temp_yz2));
dist5 = zeros(size(temp_dist5));
for i = 1:length(temp_yz2)
    if(temp_yz2(i)==prev)
        dist5(pos,:) = dist5(pos,:)+temp_dist5(i,:);
    else
        prev = temp_yz2(i);
        pos = pos+1;
        yz2(pos) = temp_yz2(i);
        dist5(pos,:) = dist5(pos,:)+temp_dist5(i,:);
    end
end
yz2 = yz2(1:pos);
dist5 = dist5(1:pos,:);

saveDistribution([name(1:dotloc-1) '_primary_area_vs_count' name(dotloc:end)],dist2,x2,t);
saveDistribution([name(1:dotloc-1) '_primary_intensity2_vs_count' name(dotloc:end)],dist5,yz2,t);
% normalized version
dist2_nbp = getNorm(dist2,prim_counts);
dist5_nbp = getNorm(dist5,prim_counts);
saveDistribution([name(1:dotloc-1) '_primary_area_vs_count_over_prim' name(dotloc:end)],dist2_nbp,x2,t);
saveDistribution([name(1:dotloc-1) '_primary_intensity2_vs_count_over_prim' name(dotloc:end)],dist5_nbp,yz2,t);

% need to get largest tumor dist
prev = -1;
pos = 0;
x3 = zeros(size(temp_x3));
dist3 = zeros(size(temp_dist3));
for i = 1:length(temp_x3)
    if(temp_x3(i)==prev)
        dist3(pos,:) = dist3(pos,:)+temp_dist3(i,:);
    else
        prev = temp_x3(i);
        pos = pos+1;
        x3(pos) = temp_x3(i);
        dist3(pos,:) = dist3(pos,:)+temp_dist3(i,:);
    end
end
x3 = x3(1:pos);
% turn 1 into zero
if(x3(1)==1)
    x3(1)=0;
end
dist3 = dist3(1:pos,:);

prev = -1;
pos = 0;
yz3 = zeros(size(temp_yz3));
dist6 = zeros(size(temp_dist6));
for i = 1:length(temp_yz3)
    if(temp_yz3(i)==prev)
        dist6(pos,:) = dist6(pos,:)+temp_dist6(i,:);
    else
        prev = temp_yz3(i);
        pos = pos+1;
        yz3(pos) = temp_yz3(i);
        dist6(pos,:) = dist6(pos,:)+temp_dist6(i,:);
    end
end
yz3 = yz3(1:pos);
dist6 = dist6(1:pos,:);

saveDistribution([name(1:dotloc-1) '_largest_area_vs_count' name(dotloc:end)],dist3,x3,t);
saveDistribution([name(1:dotloc-1) '_largest_intensity2_vs_count' name(dotloc:end)],dist6,yz3,t);
% normalized versions
dist3_nbf = getNorm(dist3,fish_counts);
dist6_nbf = getNorm(dist6,fish_counts);
saveDistribution([name(1:dotloc-1) '_largest_area_vs_count_over_fish' name(dotloc:end)],dist3_nbf,x3,t);
saveDistribution([name(1:dotloc-1) '_largest_intensity2_vs_count_over_fish' name(dotloc:end)],dist6_nbf,yz3,t);
dist3_nba = getNorm(dist3,active_fish);
dist6_nba = getNorm(dist6,active_fish);
saveDistribution([name(1:dotloc-1) '_largest_area_vs_count_over_active' name(dotloc:end)],dist3_nba,x3,t);
saveDistribution([name(1:dotloc-1) '_largest_intensity2_vs_count_over_active' name(dotloc:end)],dist6_nba,yz3,t);
end

function [norm] = getNorm(dist,factor)
norm = zeros(size(dist));
for i = 1:size(dist,2)
    if(factor(i)>0)
        norm(:,i) = dist(:,i)/factor(i);
    else
        norm(:,i) = 0;
    end
end
end

function [k_atten] = getAttenuation(A_pig)
% calculate attenuation due to pigmentation
if(A_pig > 0)
    k_atten = (1-10^(-A_pig))/A_pig/log(10);
    % don't want k_atten to be too low; the minimum
    % is based on low and high intensities in an
    % 8-bit image
    k_atten_min = (1-1/254)/log(254);
    k_atten = max(k_atten_min,k_atten);
else
    k_atten = 1;
end
end
