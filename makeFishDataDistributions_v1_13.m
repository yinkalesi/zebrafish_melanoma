%% makeFishDataDistributions_v1_13
%  Version 1.13
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 01/13/21
%% Version History
%  1.0: updating makePermutatedFishData_v1_1 for new data
%  1.1: changing min_survival and B091119 (maybe) removed due to time mismatch
%  1.2: 2019 data
%  1.3: clusterThresholdDay1 set to 10 data
%  1.4: clusterThresholdDay1 set to 10 data matching data removed
%  1.5: updated analysis (catches small shrinking tumors)
%  1.9: updating screening method to get rid of tumors that disappear on
%  day 3
%  1.10: usings rev2 files (removing day3 artifacts)
%  1.11: size restriction to > xb on individual tumor basis (using
%  getScreenedTumorDistributions_v4_0)
%  1.12 size restriction => 1<x<3e6
%  1.13: rerunning 1.10 to match fishDataTables sent to Isaac, changing to 
%  getScreenedTumorDistributions_v4_1


basedir = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\Continuous Tumor Growth Model\Heilmann-fish-image-analysis-mod\';
files_RAD = {...
    [basedir 'ct10_lt015_ht04_021720_nosat_B013017R1_VENTRAL_5e5_summary'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B021317R1_VENTRAL_5e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B070917R1_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B070917R2_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B081417R1_VENTRAL_1e5_summary'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B081417R2_VENTRAL_1e5_summary'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B081417R3_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B081417R4_VENTRAL_1e5_summary']};

files_NORAD = {...
    [basedir 'ct10_lt015_ht04_021720_nosat_B013017NR1_VENTRAL_5e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B021317NR1_VENTRAL_5e5_summary'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B072817NR1_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B072817NR2_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B072817NR3_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B072817NR4_VENTRAL_1e5_summary'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B090317NR1_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B090317NR2_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B090317NR3_VENTRAL_1e5_summary_2'],...
    [basedir 'ct10_lt015_ht04_021720_nosat_B090317NR4_VENTRAL_1e5_summary']};

files_NORAD_2019 = {
    [basedir 'ct10_lt015_ht04_021720_nosat_B090619NR_VENTRAL_1e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B091119NR_VENTRAL_1e6_summary']
    [basedir 'ct10_lt015_ht04_021720_nosat_B091919NR_VENTRAL_5e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B092319NR_VENTRAL_5e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B092719NR_VENTRAL_5e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B100119NR_VENTRAL_5e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B100919NR_VENTRAL_5e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B101319NR_MIXED_1e7_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B101719NR_DORSAL_5e6_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B102119NR_MIXED_1e7_summary_2']
    [basedir 'ct10_lt015_ht04_021720_nosat_B102519NR_MIXED_1e7_summary_2']};

files = cell(1,length(files_RAD)+length(files_NORAD)+length(files_NORAD_2019));
files(1:length(files_RAD)) = files_RAD;
files(length(files_RAD)+(1:length(files_NORAD))) = files_NORAD;
files(length(files_RAD)+length(files_NORAD)+(1:length(files_NORAD_2019))) = files_NORAD_2019;
% get times
max_tp = 20;
complete_times = zeros(1,max_tp);
% complete_times = {'1','3','5','7','9','11','13','15',...}; 
for ti = 1:max_tp
    complete_times(ti) = 2*ti-1;
end
batch_count = 29;
ntp_2017data = [13 6 6 6 9 4 5 5 5 5 8 8 8 8 8 8 8 8];
times1 = cell(1,batch_count);
for bi = 1:14
    times1{bi} = complete_times(1:ntp_2017data(bi));
end
% last for have a time point at day 4
for bi = 15:18
    times1{bi} = [complete_times(1:2) 4 complete_times(3:ntp_2017data(bi)-1)];
end
B090619_times = [0 1 3 5 7 9 11];
B091119_times = [0 2 3 4 5 8 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 58];
B091919_times = [0 1 3 5 7 9 11 13 15 17 19];
B092319_times = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 46];
B092719_times = [1 3 5 7 9 11];
B100119_times = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 38];
B100919_times = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 30];
B101319_times = [1 3 5 7 9 11 13 15 17];
B101719_times = [1 3 5 7 9 11 13 15 17 19 22];
B102119_times = [1 3 5 7 9 11 13 15 18];
B102519_times = [1 3 5 7 9 11 14];
times1{19} = (B090619_times);
times1{20} = (B091119_times);
times1{21} = (B091919_times);
times1{22} = (B092319_times);
times1{23} = (B092719_times);
times1{24} = (B100119_times);
times1{25} = (B100919_times);
times1{26} = (B101319_times);
times1{27} = (B101719_times);
times1{28} = (B102119_times);
times1{29} = (B102519_times);

% careful with fish_to_remove - there are two ordering systems (by immunity
% and by date). Make sure you are using the right one because the data is
% ordered by immunity so the input to the screening function should be by
% immunity
% fish_to_remove2017_old = {[7],[11],[1 11],[],[],[],...
%     [],[],[],[],[],[],...
%     [],[],[],[],[],[2]};
fish_to_remove2017 = {[7],[2 11],[4 11 14],[],[1 11 12],[2 5],...
    [3],[],[],[8],[],[],...
    [],[],[2 7],[],[],[2 4]};
true_ord = [1 9 2 10 3 4 11 12 13 14 5 6 7 8 15 16 17 18]; 
[~,reord] = sort(true_ord); % translate from date ordering to immunity ordering
times = times1([reord 19:29]);

% remove fish with tumors in vasculature
exclusion_list1 = [3 5; 4 3; 4 5; 4 6; 6 4; 7 9; 7 10; 11 1; 11 6; 11 7];
% removed for other reasons (speckle or erroneous data analysis)
exclusion_list2 = [2 5; 3 2; 10 1; 11 4; 11 5];
% removed because there are no tumors
exclusion_list3 = [2 7; 2 8; 7 11];
exclusion_list = [exclusion_list1; exclusion_list2; exclusion_list3];

fish_to_remove = cell(1,length(times1));
fish_to_remove(1:18) = fish_to_remove2017(reord);
for ff = 19:length(times1)
    ei = find(exclusion_list(:,1)==ff-18);
    fish_to_remove{ff} = exclusion_list(ei,2)';
end

% min_surv = 7;
for min_surv = 2 %1:7
    
    dist_dir = ['ct10_lt015_ht04_021720_newact_s' num2str(min_surv) '_011321'];
    combo_name = {'dist_RAD_zf2017.csv','dist_NORAD_zf2017.csv','dist_NORAD_zf2019.csv','dist_NORAD_all.csv'};
    gender_name = {'m','f','mf'};
    gender_id = {0, 1, [0 1]};
    if(~isfolder(dist_dir))
        mkdir(dist_dir);
    end
    
    radlist = 1:8;
    noradlist = 9:18;
    norad2019list = [19 21:29]; % skip B091119 because it does have data for t=1
    allList = [radlist noradlist norad2019list];
    ind_list = {radlist,noradlist,norad2019list,allList};
    
    % need to make gender permutations
    % 1) make file with both male and female, return which are male and which
    % are female for use in permutations
    for rnr = 1:length(combo_name)
        gi = 3;
        dotloc = find(combo_name{rnr}=='.');
        if(isempty(dotloc))
            dotloc = length(combo_name{rnr})+1;
        end
        name_gi = [combo_name{rnr}(1:dotloc-1) '_' gender_name{gi} combo_name{rnr}(dotloc:end)];
        disp(['Creating ' name_gi '...']);
        name_with_dir = [dist_dir '\' name_gi];
        start_toRemove = fish_to_remove(ind_list{rnr});
        [~,~,~,~,fdesc] = getScreenedTumorDistributions_v4_1(files(ind_list{rnr}),times(ind_list{rnr}),name_with_dir,start_toRemove,gender_id{gi},[0 inf],min_surv);
        
        % 2) use fdesc to make new toRemove list to make both male and female data
        % files separately
        % fish_disc: tumor_file,fid,gender,t1_area,tlast
        msel = find(fdesc(:,3)==0);
        mcount = length(msel);
        fsel = find(fdesc(:,3)==1);
        fcount = length(fsel);
        
        m_toRemove = start_toRemove;
        % add females to remove list
        for i = 1:fcount
            m_toRemove{fdesc(fsel(i),1)} = [m_toRemove{fdesc(fsel(i),1)} fdesc(fsel(i),2)];
        end
        
        f_toRemove = start_toRemove;
        % add females to remove list
        for i = 1:mcount
            f_toRemove{fdesc(msel(i),1)} = [f_toRemove{fdesc(msel(i),1)} fdesc(msel(i),2)];
        end
        
        gi = 1;
        name_gi = [combo_name{rnr}(1:dotloc-1) '_' gender_name{gi} combo_name{rnr}(dotloc:end)];
        disp(['Creating ' name_gi '...']);
        name_with_dir = [dist_dir '\' name_gi];
        [~,~,~,~,fdesc_m] = getScreenedTumorDistributions_v4_1(files(ind_list{rnr}),times(ind_list{rnr}),name_with_dir,m_toRemove,gender_id{gi},[0 inf],min_surv);
        desc_name_m = [dist_dir '\' combo_name{rnr}(1:dotloc-1) '_' gender_name{gi} '_desc.txt'];
        dlmwrite(desc_name_m,fdesc_m);
        
        gi = 2;
        name_gi = [combo_name{rnr}(1:dotloc-1) '_' gender_name{gi} combo_name{rnr}(dotloc:end)];
        disp(['Creating ' name_gi '...']);
        name_with_dir = [dist_dir '\' name_gi];
        [~,~,~,~,fdesc_f] = getScreenedTumorDistributions_v4_1(files(ind_list{rnr}),times(ind_list{rnr}),name_with_dir,f_toRemove,gender_id{gi},[0 inf],min_surv);
        desc_name_f = [dist_dir '\' combo_name{rnr}(1:dotloc-1) '_' gender_name{gi} '_desc.txt'];
        dlmwrite(desc_name_f,fdesc_f);
    end
end

% % generate permutations of fdesc based on mcount and fcount
% num_perms = 100;
% ord_orig = 1:mcount+fcount;
% 
% for np = 1:num_perms
%     arand = rand(1,mcount+fcount);
%     [~,ord] = sort(arand);
%     ordm = ord(1:mcount);
%     ordf = ord(mcount+1:end);
%     disp(ordm);
%     
%     % make toRemove files
%     newm_toRemove = start_toRemove;
%     for i = 1:fcount
%         newm_toRemove{fdesc(ordf(i),1)} = [newm_toRemove{fdesc(ordf(i),1)} fdesc(ordf(i),2)];
%     end
%     
%     newf_toRemove = start_toRemove;
%     for i = 1:mcount
%         newf_toRemove{fdesc(ordm(i),1)} = [newf_toRemove{fdesc(ordm(i),1)} fdesc(ordm(i),2)];
%     end
%     
%     % make distribution files
%     gi = 1;
%     fprintf('Creating permutation_%03i_%s%s ...\n',np,gender_name{gi}, combo_name{rnr}(dotloc:end));
%     name_with_dir = sprintf('%s\\permutation_%03i_%s_%s',dist_dir,np,gender_name{gi}, combo_name{rnr});
%     [~,~,~,~,fdesc_npm] = getScreenedTumorDistributions_v4_1(files(ind_list{rnr}),times(ind_list{rnr}),name_with_dir,newm_toRemove,gender_id{3},[0 inf],1);
%     desc_name1 = sprintf('%s\\permutation_%03i_%s_%s_desc.txt',dist_dir,np,gender_name{gi}, combo_name{rnr}(1:dotloc-1));
%     dlmwrite(desc_name1,fdesc_npm);
%     gi = 2;
%     disp(['Creating ' combo_name{rnr}(1:dotloc-1) '_' gender_name{gi} combo_name{rnr}(dotloc:end) '...']);
%     name_with_dir = sprintf('%s\\permutation_%03i_%s_%s',dist_dir,np,gender_name{gi}, combo_name{rnr});
%     [~,~,~,~,fdesc_npf] = getScreenedTumorDistributions_v4_1(files(ind_list{rnr}),times(ind_list{rnr}),name_with_dir,newf_toRemove,gender_id{3},[0 inf],1);
%     desc_name2 = sprintf('%s\\permutation_%03i_%s_%s_desc.txt',dist_dir,np,gender_name{gi}, combo_name{rnr}(1:dotloc-1));
%     dlmwrite(desc_name2,fdesc_npf);
% end


