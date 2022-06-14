%% plotTumorPercentiles_v1_4
%  Version 1.4
%  Author: Adeyinka Lesi
%  Date: 11/11/20
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Version History
%  1.0: from AnalyzeIndividualTrajectories_v2_4
%  1.2: appending zf2019 data; added cured condition (if fits cured
%  criteria, will fill in rest of timepoints with size 1)
%  1.3: using shifted immunity parameter...
%  1.4: version for paper (removed rad data marked cured; upon inspection, found they where likely tumor misassignment errors)
%  1.4_mod: for aiche presentation
%  1.4_mod2: making uncertainty range by individual effects of parameters
%  1.4_mod4: making plot with two shaded regions
%  1.9c: change of size limit

startt = clock();

conv = getConverter();
% want to compare intensity readings to area and area^1.5

distdir = 'ct10_lt015_ht04_021720_remerged_s2_040621/'; 
radname = 'dist_RAD_zf2017';
noradname = 'dist_NORAD_zf2017';
noradname2 = 'dist_NORAD_zf2019';

% not used yet
savedir = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\Continuous Tumor Growth Model\paper1_plots_011622\';
% version 1: size_thresh = 0, age_thresh = 0, large_tumor_thresh = 2e5
% version 2: size_thresh = 1e3, age_thresh = 0, large_tumor_thresh = 2e5
% version 3: size_thresh = 0, age_thresh = 300, large_tumor_thresh = 2e5
% version 4: size_thresh = 1e3, age_thresh = 300, large_tumor_thresh = 2e5
% version 5: size_thresh = 0, age_thresh = 0, large_tumor_thresh = 1e6
% version 6: size_thresh = 1e3, age_thresh = 0, large_tumor_thresh = 1e6
% version 7: size_thresh = 0, age_thresh = 300, large_tumor_thresh = 1e6
% version 8: size_thresh = 1e3, age_thresh = 300, large_tumor_thresh = 1e6

ver_store = [0, 0, 2e5
    1e3, 0, 2e5
    0, 300, 2e5
    1e3, 300, 2e5
    0, 0, 1e6
    1e3, 0, 1e6
    0, 300, 1e6
    1e3, 300, 1e6];

for version_num = 2%1:8

size_thresh = ver_store(version_num,1);
age_thresh = ver_store(version_num,2);
large_tumor_thresh = ver_store(version_num,3);


gender = {'m','f','mf'};
locations = cell(1,29);
for bi = 1:25
    locations{bi} = 'VENTRAL';
end
locations{26} = 'MIXED';
locations{27} = 'DORSAL';
locations{28} = 'MIXED';
locations{29} = 'MIXED';
radlist = 1:8;
noradlist = 9:18;
norad2019list = [19 21:29];

fish_batch_ages = [342   342   511   511   502   501   520   520   623   623   693   693   686   686   713   706   557   557   591   596   604   608   163   616   175   179   632   558   562];
% fish_batch_ages(1:18) = NaN;
% fish_batch_ages(19:29) = [591 596 604 608 163 616 175 179 632 558 562];
    
% need to load data
gdata = struct('int_rad',cell(1,length(gender)),...
    'area_rad',cell(1,length(gender)),...
    'label_rad',cell(1,length(gender)),...
    'desc_rad',cell(1,length(gender)),...
    'inj_loc_rad',cell(1,length(gender)),...
    'age_rad',cell(1,length(gender)),...
    'att_rad',cell(1,length(gender)),...
    'counts_rad',cell(1,length(gender)),...
    'int_norad',cell(1,length(gender)),...
    'area_norad',cell(1,length(gender)),...
    'label_norad',cell(1,length(gender)),...
    'desc_norad',cell(1,length(gender)),...
    'inj_loc_norad',cell(1,length(gender)),...
    'age_norad',cell(1,length(gender)),...
    'att_norad',cell(1,length(gender)),...
    'counts_norad',cell(1,length(gender)),...
    'int_norad2',cell(1,length(gender)),...
    'area_norad2',cell(1,length(gender)),...
    'label_norad2',cell(1,length(gender)),...
    'desc_norad2',cell(1,length(gender)),...
    'inj_loc_norad2',cell(1,length(gender)),...
    'age_norad2',cell(1,length(gender)),...
    'att_norad2',cell(1,length(gender)),...
    'counts_norad2',cell(1,length(gender)),...
    'traj_rad',cell(1,length(gender)),...
    'reltraj_rad',cell(1,length(gender)),...
    'log_rad',cell(1,length(gender)),...
    'trad',cell(1,length(gender)),...
    'rad_count',cell(1,length(gender)),...
    'traj_norad',cell(1,length(gender)),...
    'reltraj_norad',cell(1,length(gender)),...
    'log_norad',cell(1,length(gender)),...
    'tnorad',cell(1,length(gender)),...
    'norad_count',cell(1,length(gender)),...
    'traj_norad2',cell(1,length(gender)),...
    'reltraj_norad2',cell(1,length(gender)),...
    'log_norad2',cell(1,length(gender)),...
    'tnorad2',cell(1,length(gender)),...
    'norad_count2',cell(1,length(gender)),...
    'viable_rad',cell(1,length(gender)),...
    'viable_norad',cell(1,length(gender)),...
    'viable_norad2',cell(1,length(gender)),...
    'tout_r',cell(1,length(gender)),...
    'log_avg_r',cell(1,length(gender)),...
    'avg_count_r',cell(1,length(gender)),...
    'log_val_r',cell(1,length(gender)),...
    'log_med_r',cell(1,length(gender)),...
    'log_std_r',cell(1,length(gender)),...
    'tout_nr',cell(1,length(gender)),...
    'log_avg_nr',cell(1,length(gender)),...
    'avg_count_nr',cell(1,length(gender)),...
    'log_val_nr',cell(1,length(gender)),...
    'log_med_nr',cell(1,length(gender)),...   
    'log_std_nr',cell(1,length(gender)),...
    'tout_nr2',cell(1,length(gender)),...
    'log_avg_nr2',cell(1,length(gender)),...
    'avg_count_nr2',cell(1,length(gender)),...
    'log_val_nr2',cell(1,length(gender)),...
    'log_med_nr2',cell(1,length(gender)),...   
    'log_std_nr2',cell(1,length(gender)),...
    'tout_ra',cell(1,length(gender)),...
    'rel_avg_r',cell(1,length(gender)),...
    'avg_count_ra',cell(1,length(gender)),...
    'rel_val_r',cell(1,length(gender)),...
    'rel_med_r',cell(1,length(gender)),...
    'rel_std_r',cell(1,length(gender)),...
    'tout_nra',cell(1,length(gender)),...
    'rel_avg_nr',cell(1,length(gender)),...
    'avg_count_nra',cell(1,length(gender)),...
    'rel_val_nr',cell(1,length(gender)),...
    'rel_med_nr',cell(1,length(gender)),...   
    'rel_std_nr',cell(1,length(gender)),...
    'tout_nr2a',cell(1,length(gender)),...
    'rel_avg_nr2',cell(1,length(gender)),...
    'avg_count_nr2a',cell(1,length(gender)),...
    'rel_val_nr2',cell(1,length(gender)),...
    'rel_med_nr2',cell(1,length(gender)),...   
    'rel_std_nr2',cell(1,length(gender)),...
    'cured_r',cell(1,length(gender)),...
    'cured_nr',cell(1,length(gender)),...
    'cured_nr2',cell(1,length(gender)));
% male data
g = 1;
gid = gender{g};
gdata(g).int_rad = dlmread([distdir radname '_' gid '_intensities2.txt']);
gdata(g).area_rad = dlmread([distdir radname '_' gid '_sizes.txt']);
gdata(g).label_rad = dlmread([distdir radname '_' gid '_labels.txt']); % label format: f_id,tum_id,is_prim
gdata(g).desc_rad = dlmread([distdir radname '_' gid '_desc.txt']); % desc format: itf,fid,f_gender,t1_area,tlast
gdata(g).att_rad = dlmread([distdir radname '_' gid '_attenuation_values.txt']);
gdata(g).counts_rad = dlmread([distdir radname '_' gid '_counts.txt']);
gdata(g).int_norad = dlmread([distdir noradname '_' gid '_intensities2.txt']);
gdata(g).area_norad = dlmread([distdir noradname '_' gid '_sizes.txt']);
gdata(g).label_norad = dlmread([distdir noradname '_' gid '_labels.txt']);
gdata(g).desc_norad = dlmread([distdir noradname '_' gid '_desc.txt']);
gdata(g).att_norad = dlmread([distdir noradname '_' gid '_attenuation_values.txt']);
gdata(g).counts_norad = dlmread([distdir noradname '_' gid '_counts.txt']);
gdata(g).int_norad2 = dlmread([distdir noradname2 '_' gid '_intensities2.txt']);
gdata(g).area_norad2 = dlmread([distdir noradname2 '_' gid '_sizes.txt']);
gdata(g).label_norad2 = dlmread([distdir noradname2 '_' gid '_labels.txt']);
gdata(g).desc_norad2 = dlmread([distdir noradname2 '_' gid '_desc.txt']);
gdata(g).att_norad2 = dlmread([distdir noradname2 '_' gid '_attenuation_values.txt']);
gdata(g).counts_norad2 = dlmread([distdir noradname2 '_' gid '_counts.txt']);
% female data
g = 2;
gid = gender{g};
gdata(g).int_rad = dlmread([distdir radname '_' gid '_intensities2.txt']);
gdata(g).area_rad = dlmread([distdir radname '_' gid '_sizes.txt']);
gdata(g).label_rad = dlmread([distdir radname '_' gid '_labels.txt']); % label format: f_id,tum_id,is_prim
gdata(g).desc_rad = dlmread([distdir radname '_' gid '_desc.txt']); % desc format: itf,fid,f_gender,t1_area,tlast
gdata(g).att_rad = dlmread([distdir radname '_' gid '_attenuation_values.txt']);
gdata(g).counts_rad = dlmread([distdir radname '_' gid '_counts.txt']);
gdata(g).int_norad = dlmread([distdir noradname '_' gid '_intensities2.txt']);
gdata(g).area_norad = dlmread([distdir noradname '_' gid '_sizes.txt']);
gdata(g).label_norad = dlmread([distdir noradname '_' gid '_labels.txt']);
gdata(g).desc_norad = dlmread([distdir noradname '_' gid '_desc.txt']);
gdata(g).att_norad = dlmread([distdir noradname '_' gid '_attenuation_values.txt']);
gdata(g).counts_norad = dlmread([distdir noradname '_' gid '_counts.txt']);
gdata(g).int_norad2 = dlmread([distdir noradname2 '_' gid '_intensities2.txt']);
gdata(g).area_norad2 = dlmread([distdir noradname2 '_' gid '_sizes.txt']);
gdata(g).label_norad2 = dlmread([distdir noradname2 '_' gid '_labels.txt']);
gdata(g).desc_norad2 = dlmread([distdir noradname2 '_' gid '_desc.txt']);
gdata(g).att_norad2 = dlmread([distdir noradname2 '_' gid '_attenuation_values.txt']);
gdata(g).counts_norad2 = dlmread([distdir noradname2 '_' gid '_counts.txt']);
% combine data together
g = 3;
gdata(g).int_rad = [gdata(1).int_rad; gdata(2).int_rad];
gdata(g).area_rad = [gdata(1).area_rad; gdata(2).area_rad];
gdata(g).label_rad = [gdata(1).label_rad; gdata(2).label_rad]; % label format: f_id,tum_id,is_prim
gdata(g).desc_rad = [gdata(1).desc_rad; gdata(2).desc_rad]; % desc format: itf,fid,f_gender,t1_area,tlast
gdata(g).att_rad = [gdata(1).att_rad; gdata(2).att_rad];
gdata(g).counts_rad = [gdata(1).counts_rad; gdata(2).counts_rad];
gdata(g).int_norad = [gdata(1).int_norad; gdata(2).int_norad];
gdata(g).area_norad = [gdata(1).area_norad; gdata(2).area_norad];
gdata(g).label_norad = [gdata(1).label_norad; gdata(2).label_norad];
gdata(g).desc_norad = [gdata(1).desc_norad; gdata(2).desc_norad];
gdata(g).att_norad = [gdata(1).att_norad; gdata(2).att_norad];
gdata(g).counts_norad = [gdata(1).counts_norad; gdata(2).counts_norad];
gdata(g).int_norad2 = [gdata(1).int_norad2; gdata(2).int_norad2];
gdata(g).area_norad2 = [gdata(1).area_norad2; gdata(2).area_norad2];
gdata(g).label_norad2 = [gdata(1).label_norad2; gdata(2).label_norad2];
gdata(g).desc_norad2 = [gdata(1).desc_norad2; gdata(2).desc_norad2];
gdata(g).att_norad2 = [gdata(1).att_norad2; gdata(2).att_norad2];
gdata(g).counts_norad2 = [gdata(1).counts_norad2; gdata(2).counts_norad2];

tr = 1:2:25;
tnr = [1 3 4 5 7 9 11 13];
tnr2 = [1 3 5 7 9 11 13 14 15 17 18 19 21 22 23 25 27 29 30 31 33 35 37 38 39 41 43 46]; % not using day0
for g = 1:length(gender)
    gid = gender{g};
    int_rad = gdata(g).int_rad;
    area_rad = gdata(g).area_rad;
    label_rad = gdata(g).label_rad;
    desc_rad = gdata(g).desc_rad;
    att_rad = gdata(g).att_rad;
    counts_rad = gdata(g).counts_rad;
    
    int_norad = gdata(g).int_norad;
    area_norad = gdata(g).area_norad;
    label_norad = gdata(g).label_norad;
    desc_norad = gdata(g).desc_norad;
    att_norad = gdata(g).att_norad;
    counts_norad = gdata(g).counts_norad;
    
    int_norad2 = gdata(g).int_norad2;
    area_norad2 = gdata(g).area_norad2;
    label_norad2 = gdata(g).label_norad2;
    desc_norad2 = gdata(g).desc_norad2;
    att_norad2 = gdata(g).att_norad2;
    counts_norad2 = gdata(g).counts_norad2;

    size_rad = conv.a2v(area_rad);
    size_norad = conv.a2v(area_norad);
    size_norad2 = conv.a2v(area_norad2(:,2:end));
    num_sizes_rad = sum(size_rad>0,2);
    num_sizes_norad = sum(size_norad>0,2);
    num_sizes_norad2 = sum(size_norad2>0,2);
    
    % figure out location label for each tumor
    inj_loc_rad = assignLoc(locations(radlist),label_rad,desc_rad);
    inj_loc_norad = assignLoc(locations(noradlist),label_norad,desc_norad);
    inj_loc_norad2 = assignLoc(locations(norad2019list),label_norad2,desc_norad2);
    
    age_rad = assignAge(fish_batch_ages(radlist),label_rad,desc_rad);
    age_norad = assignAge(fish_batch_ages(noradlist),label_norad,desc_norad);
    age_norad2 = assignAge(fish_batch_ages(norad2019list),label_norad2,desc_norad2);
    
    % potential selection criteria: 
    % 1) at least two data points (num_sizes_rad>1) 
    % 2) is primary (label_rad(:,3)==1)
    % 3) present at time 1 (size_rad(:,1)>0)
    % 4) at least 3 data points (num_sizes_rad>2)
    viable_rad = find((num_sizes_rad>1)&(size_rad(:,1)>size_thresh));
    viable_norad = find((num_sizes_norad>1)&(size_norad(:,1)>size_thresh)); % rates vary too strongly with time/only look at tumor from Day 1
    % only interested in ventrally injected fish for zf2019
    ventral_inj_norad2 = ismember(inj_loc_norad2,{'VENTRAL'});
    gt_thresh_age_norad2 = age_norad2>age_thresh;
    viable_norad2 = find((num_sizes_norad2>1)&(size_norad2(:,1)>size_thresh)&ventral_inj_norad2&gt_thresh_age_norad2); % rates vary too strongly with time/only look at tumor from Day 1
    
    % getVal is modified getRel that gives the value, not ratio
    [traj_rad,log_rad,trad,rad_count,cured_r] = getVal(tr,viable_rad,size_rad);
    [traj_norad,log_norad,tnorad,norad_count,cured_nr] = getVal(tnr,viable_norad,size_norad);
    [traj_norad2,log_norad2,tnorad2,norad_count2,cured_nr2] = getVal(tnr2,viable_norad2,size_norad2);
    
    % Important note on method of averaging: There are 'empty data' spots 
    % in log_rad/log_norad since there isn't data at every time for every 
    % fish. The trad/tnorad variables show when there is data available
    % since trad and tnorad >=1 when there is data and are zero when there
    % is no data.
    
    % average of log
    [tout_r,log_avg_r,avg_count_r,log_val_r,log_med_r,log_std_r] = getAvg(trad(~cured_r,:),log_rad(~cured_r,:));
    [tout_nr,log_avg_nr,avg_count_nr,log_val_nr,log_med_nr,log_std_nr] = getAvg(tnorad,log_norad);
    [tout_nr2,log_avg_nr2,avg_count_nr2,log_val_nr2,log_med_nr2,log_std_nr2] = getAvg(tnorad2,log_norad2);
    
    % average of reltraj
    reltraj_rad = log10(traj_rad./traj_rad(:,1));
    reltraj_norad = log10(traj_norad./traj_norad(:,1));
    reltraj_norad2 = log10(traj_norad2./traj_norad2(:,1));
    [tout_ra,rel_avg_r,avg_count_ra,rel_val_r,rel_med_r,rel_std_r,rel_low_r,rel_high_r] = getAvg(trad(~cured_r,:),reltraj_rad(~cured_r,:));
    [tout_nra,rel_avg_nr,avg_count_nra,rel_val_nr,rel_med_nr,rel_std_nr,rel_low_nr,rel_high_nr] = getAvg(tnorad,reltraj_norad);
    [tout_nr2a,rel_avg_nr2,avg_count_nr2a,rel_val_nr2,rel_med_nr2,rel_std_nr2,rel_low_nr2,rel_high_nr2] = getAvg(tnorad2,reltraj_norad2);
    
    % try to perform t-test
    
    testres = zeros(2,length(tout_ra));
    pvals = zeros(5,length(tout_ra)); 
    for j = 2:length(tout_ra)
        tj = tout_ra(j);
        jj = find(tout_nra==tj);
        pvals(3,j) = avg_count_ra(j);
        if(~isempty(jj))
            dat1 = rel_val_r(1:avg_count_ra(j),j);
            dat2 = rel_val_nr(1:avg_count_nra(jj),jj);
            [testres(1,j),pvals(1,j)] = ttest2(dat1,dat2);
            pvals(4,j) = avg_count_nra(jj);
        end
        kk = find(tout_nr2a==tj);
        if(~isempty(kk))
            dat1 = rel_val_r(1:avg_count_ra(j),j);
            dat3 = rel_val_nr2(1:avg_count_nr2a(kk),kk);
            [testres(2,j),pvals(2,j)] = ttest2(dat1,dat3);
            pvals(5,j) = avg_count_nr2a(kk);
        end
    end

    % store results
    gdata(g).viable_rad = viable_rad;
    gdata(g).viable_norad = viable_norad;
    gdata(g).viable_norad2 = viable_norad2;
    gdata(g).inj_loc_rad = inj_loc_rad;
    gdata(g).inj_loc_norad = inj_loc_norad;
    gdata(g).inj_loc_norad2 = inj_loc_norad2;
    gdata(g).age_rad = age_rad;
    gdata(g).age_norad = age_norad;
    gdata(g).age_norad2 = age_norad2;
    
    gdata(g).traj_rad = traj_rad;
    gdata(g).reltraj_rad = reltraj_rad;
    gdata(g).log_rad = log_rad;
    gdata(g).trad = trad;
    gdata(g).rad_count = rad_count;
    gdata(g).traj_norad = traj_norad;
    gdata(g).reltraj_norad = reltraj_norad;
    gdata(g).log_norad = log_norad;
    gdata(g).tnorad = tnorad;
    gdata(g).norad_count = norad_count;
    gdata(g).traj_norad2 = traj_norad2;
    gdata(g).reltraj_norad2 = reltraj_norad2;
    gdata(g).log_norad2 = log_norad2;
    gdata(g).tnorad2 = tnorad2;
    gdata(g).norad_count2 = norad_count2;
    gdata(g).tout_r = tout_r;
    gdata(g).log_avg_r = log_avg_r;
    gdata(g).avg_count_r = avg_count_r;
    gdata(g).log_val_r = log_val_r;
    gdata(g).log_med_r = log_med_r;
    gdata(g).log_std_r = log_std_r;
    gdata(g).tout_nr = tout_nr;
    gdata(g).log_avg_nr = log_avg_nr;
    gdata(g).avg_count_nr = avg_count_nr;
    gdata(g).log_val_nr = log_val_nr;
    gdata(g).log_med_nr = log_med_nr;
    gdata(g).log_std_nr = log_std_nr;
    gdata(g).tout_nr2 = tout_nr2;
    gdata(g).log_avg_nr2 = log_avg_nr2;
    gdata(g).avg_count_nr2 = avg_count_nr2;
    gdata(g).log_val_nr2 = log_val_nr2;
    gdata(g).log_med_nr2 = log_med_nr2;
    gdata(g).log_std_nr2 = log_std_nr2;
    gdata(g).tout_ra = tout_ra;
    gdata(g).rel_avg_r = rel_avg_r;
    gdata(g).avg_count_ra = avg_count_ra;
    gdata(g).rel_val_r = rel_val_r;
    gdata(g).rel_med_r = rel_med_r;
    gdata(g).rel_std_r = rel_std_r;
    gdata(g).tout_nra = tout_nra;
    gdata(g).rel_avg_nr = rel_avg_nr;
    gdata(g).avg_count_nra = avg_count_nra;
    gdata(g).rel_val_nr = rel_val_nr;
    gdata(g).rel_med_nr = rel_med_nr;
    gdata(g).rel_std_nr = rel_std_nr;
    gdata(g).tout_nr2a = tout_nr2a;
    gdata(g).rel_avg_nr2 = rel_avg_nr2;
    gdata(g).avg_count_nr2a = avg_count_nr2a;
    gdata(g).rel_val_nr2 = rel_val_nr2;
    gdata(g).rel_med_nr2 = rel_med_nr2;
    gdata(g).rel_std_nr2 = rel_std_nr2;
    
    gdata(g).cured_r = cured_r;
    gdata(g).cured_nr = cured_nr;
    gdata(g).cured_nr2 = cured_nr2;
    
    sample_size_min = 1;
    
    % CALCULATING ERROR BARS
    % a) want to convert std of values to std of mean
    % b) we assume the relative size values have a log-normal distribuion
    
    % errors from calculating standard err from std
    devOfLogMean_low_r = rel_std_r./sqrt(avg_count_ra);
    rel_low2_r = rel_med_r-devOfLogMean_low_r;
    devOfLogMean_high_r = rel_std_r./sqrt(avg_count_ra);
    rel_high2_r = rel_med_r+devOfLogMean_high_r;
    devOfLogMean_low_nr = rel_std_nr./sqrt(avg_count_nra);
    rel_low2_nr = rel_med_nr-devOfLogMean_low_nr;
    devOfLogMean_high_nr = rel_std_nr./sqrt(avg_count_nra);
    rel_high2_nr = rel_med_nr+devOfLogMean_high_nr;
    devOfLogMean_low_nr2 = rel_std_nr2./sqrt(avg_count_nr2a);
    rel_low2_nr2 = rel_med_nr2-devOfLogMean_low_nr2;
    devOfLogMean_high_nr2 = rel_std_nr2./sqrt(avg_count_nr2a);
    rel_high2_nr2 = rel_med_nr2+devOfLogMean_high_nr2;
    
    %     % errors from high and low values from distribution
    %     devOfLogMean_low_r = (rel_med_r-rel_low_r)./sqrt(avg_count_ra);
    %     rel_low2_r = rel_med_r-devOfLogMean_low_r;
    %     devOfLogMean_high_r = (rel_high_r-rel_med_r)./sqrt(avg_count_ra);
    %     rel_high2_r = rel_med_r+devOfLogMean_high_r;
    %     devOfLogMean_low_nr = (rel_med_nr-rel_low_nr)./sqrt(avg_count_nra);
    %     rel_low2_nr = rel_med_nr-devOfLogMean_low_nr;
    %     devOfLogMean_high_nr = (rel_high_nr-rel_med_nr)./sqrt(avg_count_nra);
    %     rel_high2_nr = rel_med_nr+devOfLogMean_high_nr;
    %     devOfLogMean_low_nr2 = (rel_med_nr2-rel_low_nr2)./sqrt(avg_count_nr2a);
    %     rel_low2_nr2 = rel_med_nr2-devOfLogMean_low_nr2;
    %     devOfLogMean_high_nr2 = (rel_high_nr2-rel_med_nr2)./sqrt(avg_count_nr2a);
    %     rel_high2_nr2 = rel_med_nr2+devOfLogMean_high_nr2;
    
    lowbar_r = 10.^(rel_med_r)-10.^(rel_low2_r);
    highbar_r = 10.^(rel_high2_r)-10.^(rel_med_r);
    lowbar_nr = 10.^(rel_med_nr)-10.^(rel_low2_nr);
    highbar_nr = 10.^(rel_high2_nr)-10.^(rel_med_nr);
    lowbar_nr2 = 10.^(rel_med_nr2)-10.^(rel_low2_nr2);
    highbar_nr2 = 10.^(rel_high2_nr2)-10.^(rel_med_nr2);
    
%     plotWithErrorBars(tout_ra+1,10.^(rel_med_r),lowbar_r,highbar_r,avg_count_ra,...
%         tout_nra+1,10.^(rel_med_nr),lowbar_nr,highbar_nr,avg_count_nra,...
%         tout_nr2a+1,10.^(rel_med_nr2),lowbar_nr2,highbar_nr2,avg_count_nr2a,...
%         sample_size_min,sprintf('Rel Med, %s',gid),[-2 2]);
%     print(sprintf('%sv%i_rel_med_%s',savedir,version_num,gid),'-djpeg');
    
    % errors from calculating standard err from std
    devOfLogMean_low_r = rel_std_r./sqrt(avg_count_ra);
    rel_low2_r = rel_avg_r-devOfLogMean_low_r;
    devOfLogMean_high_r = rel_std_r./sqrt(avg_count_ra);
    rel_high2_r = rel_avg_r+devOfLogMean_high_r;
    devOfLogMean_low_nr = rel_std_nr./sqrt(avg_count_nra);
    rel_low2_nr = rel_avg_nr-devOfLogMean_low_nr;
    devOfLogMean_high_nr = rel_std_nr./sqrt(avg_count_nra);
    rel_high2_nr = rel_avg_nr+devOfLogMean_high_nr;
    devOfLogMean_low_nr2 = rel_std_nr2./sqrt(avg_count_nr2a);
    rel_low2_nr2 = rel_avg_nr2-devOfLogMean_low_nr2;
    devOfLogMean_high_nr2 = rel_std_nr2./sqrt(avg_count_nr2a);
    rel_high2_nr2 = rel_avg_nr2+devOfLogMean_high_nr2;
    
    lowbar_r = 10.^(rel_avg_r)-10.^(rel_low2_r);
    highbar_r = 10.^(rel_high2_r)-10.^(rel_avg_r);
    lowbar_nr = 10.^(rel_avg_nr)-10.^(rel_low2_nr);
    highbar_nr = 10.^(rel_high2_nr)-10.^(rel_avg_nr);
    lowbar_nr2 = 10.^(rel_avg_nr2)-10.^(rel_low2_nr2);
    highbar_nr2 = 10.^(rel_high2_nr2)-10.^(rel_avg_nr2);
%     
%     plotWithErrorBars(tout_ra+1,10.^(rel_avg_r),lowbar_r,highbar_r,avg_count_ra,...
%         tout_nra+1,10.^(rel_avg_nr),lowbar_nr,highbar_nr,avg_count_nra,...
%         tout_nr2a+1,10.^(rel_avg_nr2),lowbar_nr2,highbar_nr2,avg_count_nr2a,...
%         sample_size_min,sprintf('Rel Avg, %s',gid),[-2 2]);
%     print(sprintf('%sv%i_rel_avg_%s',savedir,version_num,gid),'-djpeg');
    
    %     dlmwrite(sprintf('%sv%i_pvals_rel_%s.txt',savedir,version_num,gid),[tout_ra+1;pvals]);

end

% plotHistogram2(gdata,0,'Male/Female ImmSupp');
% plotHistogram2(gdata,1,'Male/Female ImmComp');
% plotHistogram2(gdata,2,'Male/Female ImmComp');
% plotHistogramBySize(gdata,0,'M/F ImmSupp',[1 3 5 7 9 11 13],[1 inf]);
% plotHistogramBySize(gdata,1,'M/F ImmComp',[1 3 5 7 9 11 13],[1 inf]);
% plotHistogramBySize(gdata,2,'M/F ImmComp',[1 3 5 7 9 11 13],[1 inf]);
% plotScatter(gdata,1,'Male and Female Immunocompetent',[1:2 4:7]);
% plotScatter(gdata,0,'Male and Female Immunosuppressed',[1:7]);
% plotQuartiles(gdata,1,'Male and Female Immunocompetent',[1:2 4:7]);
% print(sprintf('%squartiles_norad',savedir),'-djpeg');
% plotQuartiles(gdata,0,'Male and Female Immunosuppressed',1:7);
% print(sprintf('%squartiles_rad',savedir),'-djpeg');
% load data with results
modresP2 = load('ct10_lt015_ht04_021720_below3e6_shiftedImmStart_sameShift_runModelHelper_zf_fits_v2_3.mat');
mkey_nr = modresP2.mrs{4,1}.key;
fkey_nr = modresP2.mrs{2,1}.key;
mkey_r = modresP2.mrs{3,1}.key;
fkey_r = modresP2.mrs{1,1}.key;
mkey_nr2 = modresP2.mrs{6,2}.key;
fkey_nr2 = modresP2.mrs{5,2}.key;

modres_hilo =load('ct10_lt015_ht04_021720_final_hilo_mm_011422_mod_runModelHelper_mm_sensitivity_v1_2.mat');
hilo_f_nr = modres_hilo.mrs(2,:);
hilo_m_nr = modres_hilo.mrs(4,:);

modres_hilo2 =load('ct10_lt015_ht04_021720_final_hilo_mm_day11Immunity_011422_mod_runModelHelper_mm_sensitivity_v1_2b.mat');
hilo_f_nr2 = modres_hilo2.mrs(5,:);
hilo_m_nr2 = modres_hilo2.mrs(6,:);
% plotWithTrajectoriesBySize2(gdata,2,'Male Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[2e5 inf],0)
% print(sprintf('%sabove_2e5_zf19_m_sep',savedir),'-djpeg');
% plotWithTrajectoriesBySize2(gdata,2,'Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[2e5 inf],1)
% print(sprintf('%sabove_2e5_zf19_f_sep',savedir),'-djpeg');

% plotWithTrajectoriesBySize3(gdata,2,'Male Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[2e5 inf],0)
% print(sprintf('%sabove_2e5_zf19_m_sep_wtheory',savedir),'-djpeg');
% plotWithTrajectoriesBySize3(gdata,2,'Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[2e5 inf],1)
% print(sprintf('%sabove_2e5_zf19_f_sep_wtheory',savedir),'-djpeg');
% 
% plotWithTrajectoriesBySize3(gdata,2,'Male Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr2,fkey_nr2,[2e5 inf],0)
% print(sprintf('%sabove_2e5_zf19_m_sep_wtheory2',savedir),'-djpeg');
% plotWithTrajectoriesBySize3(gdata,2,'Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr2,fkey_nr2,[2e5 inf],1)
% print(sprintf('%sabove_2e5_zf19_f_sep_wtheory2',savedir),'-djpeg');

% plotWithTrajectoriesBySize5(gdata,2,'',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,mkey_nr2,fkey_nr2,[large_tumor_thresh inf],0)
% print(sprintf('%sv%i_above_2e5_zf19_m_sep_wtheory3',savedir,version_num),'-djpeg');
% plotWithTrajectoriesBySize5(gdata,2,'',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,mkey_nr2,fkey_nr2,[large_tumor_thresh inf],1)
% print(sprintf('%sv%i_above_2e5_zf19_f_sep_wtheory3',savedir,version_num),'-djpeg');

% plotRelativeSize(gdata,2,'',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,mkey_nr2,fkey_nr2,[1 inf],0);
% plotRelativeSize(gdata,2,'',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,mkey_nr2,fkey_nr2,[1 inf],1);

% plotWithTrajectoriesAreaError2(gdata,2,'',[2:8 10:11 13:14]-1,hilo_m_nr,[large_tumor_thresh inf],0)
% print(sprintf('%sv%i_above_2e5_zf19_m_areaerror',savedir,version_num),'-djpeg');
% plotWithTrajectoriesAreaError2(gdata,2,'',[2:8 10:11 13:14]-1,hilo_f_nr,[large_tumor_thresh inf],1)
% print(sprintf('%sv%i_above_2e5_zf19_f_areaerror',savedir,version_num),'-djpeg');
% 
% plotWithTrajectoriesAreaError2(gdata,2,'',[2:8 10:11 13:14]-1,hilo_m_nr2,[large_tumor_thresh inf],0)
% print(sprintf('%sv%i_above_2e5_zf19_m_areaerror_shiftedImmStart',savedir,version_num),'-djpeg');
% plotWithTrajectoriesAreaError2(gdata,2,'',[2:8 10:11 13:14]-1,hilo_f_nr2,[large_tumor_thresh inf],1)
% print(sprintf('%sv%i_above_2e5_zf19_f_areaerror_shiftedImmStart',savedir,version_num),'-djpeg');

% plotWithTrajectoriesBySize(gdata,1,'',[1:2 4:7],mkey_nr,fkey_nr,[1 inf]);
% print(sprintf('%sv%i_mf_prediction_norad',savedir,version_num),'-djpeg');

plotWithTrajectoriesAreaError3(gdata,2,'',[2:8 10:11 13:14]-1,hilo_m_nr,hilo_m_nr2,[large_tumor_thresh inf],0)
print(sprintf('%sv%i_above_2e5_zf19_m_areaerror_2predictions',savedir,version_num),'-djpeg');
plotWithTrajectoriesAreaError3(gdata,2,'',[2:8 10:11 13:14]-1,hilo_f_nr,hilo_f_nr2,[large_tumor_thresh inf],1)
print(sprintf('%sv%i_above_2e5_zf19_f_areaerror_2predictions',savedir,version_num),'-djpeg');

% % day 1 size distributions
% plotHistogramDay1(gdata,0,1,'Male Immune-Suppressed');
% print(sprintf('%sv%i_hist_r_m',savedir,version_num),'-djpeg');
% plotHistogramDay1(gdata,0,2,'Female Immune-Suppressed');
% print(sprintf('%sv%i_hist_r_f',savedir,version_num),'-djpeg');
% plotHistogramDay1(gdata,1,1,'Male Immune-Competent (Batch 1)');
% print(sprintf('%sv%i_hist_nr_m',savedir,version_num),'-djpeg');
% plotHistogramDay1(gdata,1,2,'Female Immune-Competent (Batch 1)');
% print(sprintf('%sv%i_hist_nr_f',savedir,version_num),'-djpeg');
% plotHistogramDay1(gdata,2,1,'Male Immune-Competent (Batch 2)');
% print(sprintf('%sv%i_hist_nr2_m',savedir,version_num),'-djpeg');
% plotHistogramDay1(gdata,2,2,'Female Immune-Competent (Batch 2)');
% print(sprintf('%sv%i_hist_nr2_f',savedir,version_num),'-djpeg');

% % gender gender difference pvals
% g_pvals1 = zeros(3,length(gdata(3).tout_ra));
% g_test1 = zeros(1,length(gdata(3).tout_ra));
% for j = 1:length(g_test1)
%     dat1 = log10(gdata(1).rel_val_r(1:gdata(1).avg_count_ra(j),j));
%     dat2 = log10(gdata(2).rel_val_r(1:gdata(2).avg_count_ra(j),j));
%     [g_test1(1,j),g_pvals1(1,j)] = ttest2(dat1,dat2);
%     g_pvals1(2,j) = gdata(1).avg_count_ra(j);
%     g_pvals1(3,j) = gdata(2).avg_count_ra(j);
% end
% dlmwrite(sprintf('%sv%i_gender_pvals_rel_r.txt',savedir,version_num),[gdata(3).tout_ra+1;g_pvals1]);
% 
% g_pvals2 = zeros(3,length(gdata(3).tout_nra));
% g_test2 = zeros(1,length(gdata(3).tout_nra));
% for j = 1:length(g_test2)
%     dat1 = gdata(1).rel_val_nr(1:gdata(1).avg_count_nra(j),j);
%     dat2 = gdata(2).rel_val_nr(1:gdata(2).avg_count_nra(j),j);
%     [g_test2(1,j),g_pvals2(1,j)] = ttest2(dat1,dat2);
%     g_pvals2(2,j) = gdata(1).avg_count_nra(j);
%     g_pvals2(3,j) = gdata(2).avg_count_nra(j);
% end
% dlmwrite(sprintf('%sv%i_gender_pvals_nr.txt',savedir,version_num),[gdata(3).tout_nra+1;g_pvals2]);
% 
% g_pvals3 = zeros(3,length(gdata(3).tout_nr2a));
% g_test3 = zeros(1,length(gdata(3).tout_nr2a));
% for j = 1:length(g_test3)
%     dat1 = gdata(1).rel_val_nr2(1:gdata(1).avg_count_nr2a(j),j);
%     dat2 = gdata(2).rel_val_nr2(1:gdata(2).avg_count_nr2a(j),j);
%     [g_test3(1,j),g_pvals3(1,j)] = ttest2(dat1,dat2);
%     g_pvals3(2,j) = gdata(1).avg_count_nr2a(j);
%     g_pvals3(3,j) = gdata(2).avg_count_nr2a(j);
% end
% dlmwrite(sprintf('%sv%i_gender_pvals_nr2.txt',savedir,version_num),[gdata(3).tout_nr2a+1;g_pvals3]);

% plotWithTrajectoriesBySize1(gdata,1,'Male and Female Immunocompetent, 2017 Data',[1:2 4:7],mkey_nr,fkey_nr,[2e5 inf])
% print(sprintf('%sabove_2e5_norad',savedir),'-djpeg');
% % plotWithTrajectories(gdata,0,'Male and Female Immunosuppressed',1:7,mkey_r,fkey_r)
% % print(sprintf('%sabove_crossover_rad',savedir),'-djpeg');
% plotWithTrajectoriesBySize1(gdata,2,'Male and Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[2e5 inf])
% print(sprintf('%sabove_2e5_norad2',savedir),'-djpeg');

% plotWithTrajectoriesBySize1(gdata,1,'Male and Female Immunocompetent, 2017 Data',[1:2 4:7],mkey_nr,fkey_nr,[3e5 2e5])
% print(sprintf('%sabove_3e5_norad',savedir),'-djpeg');
% % plotWithTrajectories(gdata,0,'Male and Female Immunosuppressed',1:7,mkey_r,fkey_r)
% % print(sprintf('%sabove_crossover_rad',savedir),'-djpeg');
% plotWithTrajectoriesBySize1(gdata,2,'Male and Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[3e5 2e5])
% print(sprintf('%sabove_3e5_norad2',savedir),'-djpeg');
% 
% plotWithTrajectoriesBySize1(gdata,1,'Male and Female Immunocompetent, 2017 Data',[1:2 4:7],mkey_nr,fkey_nr,[3e4 3e5])
% print(sprintf('%sabove_3e4_norad',savedir),'-djpeg');
% % plotWithTrajectories(gdata,0,'Male and Female Immunosuppressed',1:7,mkey_r,fkey_r)
% % print(sprintf('%sabove_crossover_rad',savedir),'-djpeg');
% plotWithTrajectoriesBySize1(gdata,2,'Male and Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[3e4 3e5])
% print(sprintf('%sabove_3e4_norad2',savedir),'-djpeg');
% 
% plotWithTrajectoriesBySize1(gdata,1,'Male and Female Immunocompetent, 2017 Data',[1:2 4:7],mkey_nr,fkey_nr,[3e3 3e4])
% print(sprintf('%sabove_3e3_norad',savedir),'-djpeg');
% % plotWithTrajectories(gdata,0,'Male and Female Immunosuppressed',1:7,mkey_r,fkey_r)
% % print(sprintf('%sabove_crossover_rad',savedir),'-djpeg');
% plotWithTrajectoriesBySize1(gdata,2,'Male and Female Immunocompetent, 2019 Data',[2:8 10:11 13:14]-1,mkey_nr,fkey_nr,[3e3 3e4])
% print(sprintf('%sabove_3e3_norad2',savedir),'-djpeg');

    
save(sprintf('%sv%i_gdata.mat',savedir,version_num),'gdata');
    
close all;
end

endt = clock();
duramins = etime(endt,startt)/60;

fprintf('%s took %3.2f minutes\n',mfilename,duramins);

function [reltraj,log_size,tout,count] = getRel(t,viable,sizes)
reltraj = zeros(length(viable),length(t));
tout = zeros(length(viable),length(t));
log_size = zeros(length(viable),length(t));
count = zeros(1,length(t));

for i = 1:length(viable)
    ii = viable(i);
    traj = sizes(ii,:);
    js = find(traj~=0);
    if(~isempty(js))
        reltraj(i,1:length(js)) = traj(js)/traj(js(1));
        tout(i,1:length(js)) = t(js);
%         plot(tout(i,1:length(js)),reltraj(i,1:length(js)));
        log_size(i,1:length(js)) = log(reltraj(i,1:length(js)));
        count(js) = count(js)+1;
%         hold on;
    end
end
% hold off;
end

function [trajs,log_size,tout,count,cured] = getVal(t,viable,sizes)
trajs = zeros(length(viable),length(t));
tout = zeros(length(viable),length(t));
log_size = zeros(length(viable),length(t));
count = zeros(1,length(t));
cured = zeros(length(viable),1);
CURED_THRES = 20;
LOW_SIZE_FILL = 4.0167;

for i = 1:length(viable)
    ii = viable(i);
    traj = sizes(ii,:);
    js = find(traj~=0);
    if(~isempty(js))
        cured(i) = traj(js(end)) <= CURED_THRES;
        count(js) = count(js)+1;
        if(cured(i))
            % set rest of time points to LOW_SIZE_FILL if cured
            trajs(i,1:length(js)) = traj(js);
            trajs(i,length(js)+(1:length(t)-js(end))) = LOW_SIZE_FILL;
            tout(i,1:length(js)) = t(js);
            tout(i,length(js)+(1:length(t)-js(end))) = t(js(end)+1:end);
        else
            trajs(i,1:length(js)) = traj(js);
            tout(i,1:length(js)) = t(js);
        end
%         plot(tout(i,1:length(js)),reltraj(i,1:length(js)));
        log_size(i,trajs(i,:)~=0) = log(trajs(i,trajs(i,:)~=0));
%         hold on;
    end
end
% hold off;
end

function [] = plotRel(ts,rels)
for i = 1:size(rels,1)
    nz = find(rels(i,:)~=0);
    semilogy(ts(i,nz)-ts(i,nz(1)),rels(i,nz));
    hold on;
end
hold off;
end

function [tout,log_ave,count,log_val,log_med,log_std,low34,high34] = getAvg(ts,log_size)
ts2 = ts-ts(:,1);
tout = unique(ts2)';
tout = tout(find(tout==0,1):end);
log_ave = zeros(1,length(tout));
log_val = zeros(size(log_size,1),length(tout));
log_std = zeros(1,length(tout));
low34 = zeros(1,length(tout));
high34 = zeros(1,length(tout));
log_med = zeros(1,length(tout));
count = zeros(1,length(tout));
lookup = zeros(1,tout(end)-tout(1)+1);
for k = 1:length(tout)
    lookup(tout(k)-tout(1)+1) = k;
end
for i = 1:size(ts2,1)
    nn = find(ts2(i,:)>=0);
    ind = lookup(ts2(i,nn)-tout(1)+1);
    log_ave(ind) = log_ave(ind)+log_size(i,nn);
    count(ind) = count(ind)+1;
    for j = 1:length(ind)
        log_val(count(ind(j)),ind(j)) = log_size(i,nn(j));
    end
end
nz = find(count~=0);
log_ave(nz) = log_ave(nz)./count(nz);
for i = 1:length(nz)
    ii = nz(i);
    if(ii > 0) % first time point is included for this script
        vals_ii = sort(log_val(1:count(ii),ii));
        log_std(ii) = std(vals_ii);
        log_med(ii) = median(vals_ii);
        % calculate values 34% below median and 34% above as an alternative
        % to using standard deviation
        low34_ind = max(1,round(count(ii)*(0.5-0.34)));
        low34(ii) = vals_ii(low34_ind);
        high34_ind = round(count(ii)*(0.5+0.34));
        high34(ii) = vals_ii(high34_ind);
    end
end
end

function [] = plotWithErrorBars(tout_r,avg_r,lowbar_r,highbar_r,avg_count_r,...
    tout_nr,avg_nr,lowbar_nr,highbar_nr,avg_count_nr,...
    tout_nr2,avg_nr2,lowbar_nr2,highbar_nr2,avg_count_nr2,...
    sample_size_min,title_text,log_ylims)
legpos1 = [0.7 0.55 0.2 0.4];
legpos2 = [0.75 0.45 0.2 0.4];
legpos3 = [0.75 0.45 0.2 0.4];
xtic1 = 10.^(7:1:12);
xtic2 = 10.^(0:1:5);
xtic3 = 10.^(0:1:9);
ytic_def = 'auto';
hpic = 3.125*61*2;
wpic = 3.125*61*2;
fontsizes = [6 8 7 6]*2; 

fig = figure;
fig.Position = [962 42 wpic hpic];
sel = find(avg_count_r>=sample_size_min);
errorbar(tout_r(sel),avg_r(sel),lowbar_r(sel),highbar_r(sel),'o-','MarkerSize',3,'LineWidth',1.5);
hold on;
% don't want to plot data for 1 day difference (use only data with n>=limit)
sel = find(avg_count_nr>=sample_size_min);
errorbar(tout_nr(sel),avg_nr(sel),lowbar_nr(sel),highbar_nr(sel),'o-','MarkerSize',3,'LineWidth',1.5);
% second set of nr data
sel = find(avg_count_nr2>=sample_size_min);
errorbar(tout_nr2(sel),avg_nr2(sel),lowbar_nr2(sel),highbar_nr2(sel),'o-','MarkerSize',3,'LineWidth',1.5);
hold off;
set(gca,'YScale','log');
axis([1 10 10.^(log_ylims)]);
yticks(10.^(log_ylims(1):1:log_ylims(2)));
xlabel('Days After Inoculation','FontSize',fontsizes(2));
ylabel('Tumor Volume/Day 1 Volume','FontSize',fontsizes(3));
legend({'Immunosuppressed (Set 1)','Immunocompetent (Set 2)','Immunocompetent (Set 3)'},'Location','Best','FontSize',fontsizes(4));
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
% title(title_text);
end

function [] = plotHistogram(tout,log_val,avg_count,log_avg,log_med,desc)
figure('units','normalized','outerposition',[.05 .05 .90 .90]);
nplot = length(tout)-1; % skip first time point
ncol = ceil(sqrt(nplot));
nrow = ceil(nplot/ncol);
for i = 1:nplot
    subplot(nrow,ncol,i);
    histogram(log_val(1:avg_count(i+1),i+1),20);
    xlabel('Log Size Ratio');
    hold on;
    plot([1 1]*log_avg(i+1),[0 5],'o-',...
        [1 1]*log_med(i+1),[0 5],'s-','LineWidth',1.5);
    legend('Count','Mean','Median','Location','Best');
    title(sprintf('%s, \\Delta{t}=%i, n=%i',desc,tout(i+1),avg_count(i+1)));
    hold off;
end
end

function [] = plotHistogramDay1(gdata,rnr,gid,title_text)
figure('units','normalized','outerposition',[.05 .05 .45 .55]);

if(rnr==0)
    log_val = gdata(gid).log_rad(:,1)/log(10);
elseif(rnr==1)
    log_val = gdata(gid).log_norad(:,1)/log(10);
else
    log_val = gdata(gid).log_norad2(:,1)/log(10);
end

histogram(log_val,'BinEdges',linspace(0,7,8));
xlabel('Log Tumor Size');
title(title_text);
end

function [] = plotHistogram2(gdata,rnr,desc)
% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .90 .90]);

if(rnr == 0)
    tout1 = gdata(1).tout_r;
    log_val1 = gdata(1).log_val_r;
    avg_count1 = gdata(1).avg_count_r;
    log_avg1 = gdata(1).log_avg_r;
    log_med1 = gdata(1).log_med_r;
    tout2 = gdata(2).tout_r;
    log_val2 = gdata(2).log_val_r;
    avg_count2 = gdata(2).avg_count_r;
    log_avg2 = gdata(2).log_avg_r;
    log_med2 = gdata(2).log_med_r;
elseif(rnr == 1)
    tout1 = gdata(1).tout_nr;
    log_val1 = gdata(1).log_val_nr;
    avg_count1 = gdata(1).avg_count_nr;
    log_avg1 = gdata(1).log_avg_nr;
    log_med1 = gdata(1).log_med_nr;
    tout2 = gdata(2).tout_nr;
    log_val2 = gdata(2).log_val_nr;
    avg_count2 = gdata(2).avg_count_nr;
    log_avg2 = gdata(2).log_avg_nr;
    log_med2 = gdata(2).log_med_nr;
else
    tout1 = gdata(1).tout_nr2(1:10); % only want days with sufficient data
    log_val1 = gdata(1).log_val_nr2;
    avg_count1 = gdata(1).avg_count_nr2;
    log_avg1 = gdata(1).log_avg_nr2;
    log_med1 = gdata(1).log_med_nr2;
    tout2 = gdata(2).tout_nr2(1:10);
    log_val2 = gdata(2).log_val_nr2;
    avg_count2 = gdata(2).avg_count_nr2;
    log_avg2 = gdata(2).log_avg_nr2;
    log_med2 = gdata(2).log_med_nr2;
end
nplot = length(tout1); % skip first time point
ncol = ceil(sqrt(nplot));
nrow = ceil(nplot/ncol);
% set up histogram bins
vmin = floor(min([log_val1(:); log_val1(:)]));
vmax = ceil(max([log_val1(:); log_val1(:)]));
edges = linspace(vmin,vmax,10);
for i = 1:nplot
    subplot(nrow,ncol,i);
    
    val1 = log_val1(1:avg_count1(i),i);
    val2 = log_val2(1:avg_count2(i),i);
    count1 = histcounts(val1,edges);
    count2 = histcounts(val2,edges);
    dist1 = count1/length(val1);
    dist2 = count2/length(val2);
    mids = 0.5*(edges(2:end)+edges(1:end-1));
    bar(mids,[dist1' dist2'],'BarWidth',1.5);
    xlabel('Log Size');
%     legend('M Count','F Count','Location','Best');
    hold on;
    plot([1 1]*log_avg1(i),[0 0.2],'o-',[1 1]*log_med1(i),[0 0.2],'s-',...
         [1 1]*log_avg2(i),[0 0.2],'o-',[1 1]*log_med2(i),[0 0.2],'s-','LineWidth',1.0);
    legend('M Count','F Count','M Mean','M Median','F Mean','F Median','Location','Best');
    title(sprintf('%s, \\Delta{t}=%i, m=%i, f=%i',desc,tout1(i),avg_count1(i),avg_count2(i)));
    hold off;
end
end

function [] = plotHistogramBySize(gdata,rnr,desc,tplot,size_range)
% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .90 .90]);

if(rnr == 0)
    tout1 = gdata(1).tout_r;
    t1 = gdata(1).trad;
    log_val1 = gdata(1).log_rad;
    avg_count1 = gdata(1).avg_count_r;
    log_avg1 = gdata(1).log_avg_r;
    log_med1 = gdata(1).log_med_r;
    tout2 = gdata(2).tout_r;
    t2 = gdata(2).trad;
    log_val2 = gdata(2).log_rad;
    avg_count2 = gdata(2).avg_count_r;
    log_avg2 = gdata(2).log_avg_r;
    log_med2 = gdata(2).log_med_r;
elseif(rnr == 1)
    tout1 = gdata(1).tout_nr;
    t1 = gdata(1).tnorad;
    log_val1 = gdata(1).log_norad;
    avg_count1 = gdata(1).avg_count_nr;
    log_avg1 = gdata(1).log_avg_nr;
    log_med1 = gdata(1).log_med_nr;
    tout2 = gdata(2).tout_nr;
    t2 = gdata(2).tnorad;
    log_val2 = gdata(2).log_norad;
    avg_count2 = gdata(2).avg_count_nr;
    log_avg2 = gdata(2).log_avg_nr;
    log_med2 = gdata(2).log_med_nr;
else
    tout1 = gdata(1).tout_nr2;
    t1 = gdata(1).tnorad2;
    log_val1 = gdata(1).log_norad2;
    avg_count1 = gdata(1).avg_count_nr2;
    log_avg1 = gdata(1).log_avg_nr2;
    log_med1 = gdata(1).log_med_nr2;
    tout2 = gdata(2).tout_nr2;
    t2 = gdata(2).tnorad2;
    log_val2 = gdata(2).log_norad2;
    avg_count2 = gdata(2).avg_count_nr2;
    log_avg2 = gdata(2).log_avg_nr2;
    log_med2 = gdata(2).log_med_nr2;
end

% select tumors within size range
in_range1 = (t1(:,1)==1)&(log_val1(:,1)>=log(size_range(1)))&(log_val1(:,1)<=log(size_range(2)));
in_range2 = (t2(:,1)==1)&(log_val2(:,1)>=log(size_range(1)))&(log_val2(:,1)<=log(size_range(2)));

in_val1 = log_val1(in_range1,:);
in_val2 = log_val2(in_range2,:);
in_t1 = t1(in_range1,:);
in_t2 = t2(in_range2,:);

nplot = length(tplot);
ncol = ceil(sqrt(nplot));
nrow = ceil(nplot/ncol);
% set up histogram bins
vmin = log(1e2);%floor(min(nonzeros([in_val1(:); in_val2(:)])));
vmax = ceil(max(nonzeros([in_val1(:); in_val2(:)])));
edges = linspace(vmin,vmax,10);
for i = 1:nplot
    subplot(nrow,ncol,i);
    
    % find location of proper quantity for each row
    ind1 = ismember(in_t1,tplot(i));
    ind2 = ismember(in_t2,tplot(i));
    ii = find(tout1+1==tplot(i));
    
    val1 = in_val1(ind1);
    val2 = in_val2(ind2);
    count1 = histcounts(val1,edges);
    count2 = histcounts(val2,edges);
    dist1 = count1/length(val1);
    dist2 = count2/length(val2);
    mids = 0.5*(edges(2:end)+edges(1:end-1));
    plot(mids,[dist1' dist2'],'o-','LineWidth',1.5);
    xlabel('Log Size');
%     legend('M Count','F Count','Location','Best');
    hold on;
    legend('M Count','F Count','Location','Best');
    %     plot([1 1]*log_avg1(ii),[0 0.2],'o-',[1 1]*log_med1(ii),[0 0.2],'s-',...
    %          [1 1]*log_avg2(ii),[0 0.2],'o-',[1 1]*log_med2(ii),[0 0.2],'s-','LineWidth',1.0);
    %     legend('M Count','F Count','M Mean','M Median','F Mean','F Median','Location','Best');
    title(sprintf('%s, \\Delta{t}=%i, m=%i, f=%i',desc,tout1(ii),avg_count1(ii),avg_count2(ii)));
    hold off;
end
end

function [] = plotScatter(gdata,rnr,desc,tsel)
% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .90 .90]);
fontsizes = [6 8 7 6]*2; 
if(rnr == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
else
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
end

% let's scatter plot and quartile
split = [0.25 0.5 0.75]';
quart1 = zeros(length(split),length(tout1));
quart2 = zeros(length(split),length(tout2));
for i = 1:length(tout1)
    mh1 = plot((tout1(i)+1-.05)*ones(1,avg_count1(i)),log_val1(1:avg_count1(i),i)/log(10),'bx');
    hold on;
    mh2 = plot((tout2(i)+1.05)*ones(1,avg_count2(i)),log_val2(1:avg_count2(i),i)/log(10),'r+');
    % get quartiles
    [sval1,~] = sort(log_val1(1:avg_count1(i),i));
    [sval2,~] = sort(log_val2(1:avg_count2(i),i));
    iq1_lo = 1+floor(split*(avg_count1(i)-1));
    iq1_hi = 1+ceil(split*(avg_count1(i)-1));
    iq2_lo = 1+floor(split*(avg_count2(i)-1));
    iq2_hi = 1+ceil(split*(avg_count2(i)-1));
    nq1_lo = (iq1_lo-1)/(avg_count1(i)-1);
    nq1_hi = (iq1_hi-1)/(avg_count1(i)-1);
    nq2_lo = (iq2_lo-1)/(avg_count2(i)-1);
    nq2_hi = (iq2_hi-1)/(avg_count2(i)-1);
    slp1 = (sval1(iq1_hi)-sval1(iq1_lo))./(nq1_hi-nq1_lo);
    slp1(nq1_hi-nq1_lo==0) = 0;
    quart1(:,i) = sval1(iq1_lo)+slp1.*(split-nq1_lo);
    slp2 = (sval2(iq2_hi)-sval2(iq2_lo))./(nq2_hi-nq2_lo);
    slp2(nq2_hi-nq2_lo==0) = 0;
    quart2(:,i) = sval2(iq2_lo)+slp2.*(split-nq2_lo);
end

h = plot(tout1+st,quart1/log(10),'b-',tout2+st,quart2/log(10),'r-',...
    'LineWidth',1.5);
set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
% legend('Male','Female');
legend([mh1; mh2; h],'Male, Data','Female, Data','Male, 1st Quartile',...
    'Male, Median','Male, 3rd Quartile','Female, 1st Quartile',...
    'Female, Median','Female, 3rd Quartile','Location','Best');
xlabel('Days after Inoculation','FontSize',fontsizes(2));
ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
title(sprintf('%s',desc));
% get detection limit;
% conv = getConverter();
% dlim = conv.a2v(3);
dlim = 100;
axis([st-0.2 st+tout1(end)+0.2 log10(dlim) inf]);
end

function [] = plotQuartiles(gdata,rnr,desc,tsel)
% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .90 .90]);
fontsizes = [6 8 7 6]*2; 
if(rnr == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
else
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
end

% let's scatter plot and quartile
split = [0.25 0.5 0.75]';
quart1 = zeros(length(split),length(tout1));
quart2 = zeros(length(split),length(tout2));
for i = 1:length(tout1)
%     mh1 = plot((tout1(i)+1-.05)*ones(1,avg_count1(i)),log_val1(1:avg_count1(i),i)/log(10),'bx');
%     hold on;
%     mh2 = plot((tout2(i)+1.05)*ones(1,avg_count2(i)),log_val2(1:avg_count2(i),i)/log(10),'r+');
    % get quartiles
    [sval1,~] = sort(log_val1(1:avg_count1(i),i));
    [sval2,~] = sort(log_val2(1:avg_count2(i),i));
    iq1_lo = 1+floor(split*(avg_count1(i)-1));
    iq1_hi = 1+ceil(split*(avg_count1(i)-1));
    iq2_lo = 1+floor(split*(avg_count2(i)-1));
    iq2_hi = 1+ceil(split*(avg_count2(i)-1));
    nq1_lo = (iq1_lo-1)/(avg_count1(i)-1);
    nq1_hi = (iq1_hi-1)/(avg_count1(i)-1);
    nq2_lo = (iq2_lo-1)/(avg_count2(i)-1);
    nq2_hi = (iq2_hi-1)/(avg_count2(i)-1);
    slp1 = (sval1(iq1_hi)-sval1(iq1_lo))./(nq1_hi-nq1_lo);
    slp1(nq1_hi-nq1_lo==0) = 0;
    quart1(:,i) = sval1(iq1_lo)+slp1.*(split-nq1_lo);
    slp2 = (sval2(iq2_hi)-sval2(iq2_lo))./(nq2_hi-nq2_lo);
    slp2(nq2_hi-nq2_lo==0) = 0;
    quart2(:,i) = sval2(iq2_lo)+slp2.*(split-nq2_lo);
end

h = plot(tout1+st,quart1/log(10),'b-',tout2+st,quart2/log(10),'r-',...
    'LineWidth',1.5);
set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
legend(h,'Male, 1st Quartile',...
    'Male, Median','Male, 3rd Quartile','Female, 1st Quartile',...
    'Female, Median','Female, 3rd Quartile','Location','Best');
xlabel('Days after Inoculation','FontSize',fontsizes(2));
ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
title(sprintf('%s',desc));
% get detection limit;
% conv = getConverter();
% dlim = conv.a2v(3);
dlim = 100;
axlims = [st-0.5 st+tout1(end)+0.5 log10(dlim) 8];
axis(axlims);

edges = linspace(axlims(3),axlims(4),51);
ax1 = gca;
axslp = ax1.Position(3)/(axlims(2)-axlims(1));
for i = 1:length(tout1)
    pwid = 0.025;
    ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
    ax2 = axes('Position',ax2pos);
    sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
    barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
    axis([-10 0 axlims(3:4)]);
    ax2.XTickLabel = [];
    ax2.XTick = [];
    ax2.YTickLabel = [];
    ax2.YTick = [];
    ax2.Visible = 'off';
    ax2.XColor = 'none';
    ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
    ax3 = axes('Position',ax3pos);
    sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
    barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
    axis([0 10 axlims(3:4)]);
    ax3.XTickLabel = [];
    ax3.XTick = [];
    ax3.YTickLabel = [];
    ax3.YTick = [];
    ax3.Visible = 'off';
    ax3.XColor = 'none';
end
end

% incorporate model fit
function [] = plotWithTrajectories(gdata,rnr,desc,tsel,mkey,fkey)
% get size trajectory of male and females
dt = 0.01;
if rnr == 0
    tf = 13;
    x0 = 3.5e5;
    x1 = 3.5e6;
else
    tf = 12;
    x0 = 6e5;
    x1 = 6e6;
end
[mtraj0,mt0] = getSizeTrajectory_v3_3(mkey,x0,dt,tf);
lmtraj0 = log10(mtraj0);
[ftraj0,ft0] = getSizeTrajectory_v3_3(fkey,x0,dt,tf);
lftraj0 = log10(ftraj0);
[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey,x1,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey,x1,dt,tf);
lftraj1 = log10(ftraj1);
% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .75 .85]);
fontsizes = [6 8 8 7]*2; 
if(rnr == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
elseif(rnr == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
end

% let's scatter plot and quartile
split = [0.25 0.5 0.75]';
quart1 = zeros(length(split),length(tout1));
quart2 = zeros(length(split),length(tout2));
for i = 1:length(tout1)
%     mh1 = plot((tout1(i)+1-.05)*ones(1,avg_count1(i)),log_val1(1:avg_count1(i),i)/log(10),'bx');
%     hold on;
%     mh2 = plot((tout2(i)+1.05)*ones(1,avg_count2(i)),log_val2(1:avg_count2(i),i)/log(10),'r+');
    % get quartiles
    [sval1,~] = sort(log_val1(1:avg_count1(i),i));
    [sval2,~] = sort(log_val2(1:avg_count2(i),i));
    iq1_lo = 1+floor(split*(avg_count1(i)-1));
    iq1_hi = 1+ceil(split*(avg_count1(i)-1));
    iq2_lo = 1+floor(split*(avg_count2(i)-1));
    iq2_hi = 1+ceil(split*(avg_count2(i)-1));
    nq1_lo = (iq1_lo-1)/(avg_count1(i)-1);
    nq1_hi = (iq1_hi-1)/(avg_count1(i)-1);
    nq2_lo = (iq2_lo-1)/(avg_count2(i)-1);
    nq2_hi = (iq2_hi-1)/(avg_count2(i)-1);
    slp1 = (sval1(iq1_hi)-sval1(iq1_lo))./(nq1_hi-nq1_lo);
    slp1(nq1_hi-nq1_lo==0) = 0;
    quart1(:,i) = sval1(iq1_lo)+slp1.*(split-nq1_lo);
    slp2 = (sval2(iq2_hi)-sval2(iq2_lo))./(nq2_hi-nq2_lo);
    slp2(nq2_hi-nq2_lo==0) = 0;
    quart2(:,i) = sval2(iq2_lo)+slp2.*(split-nq2_lo);
end

h = plot(tout1+st,quart1(3,:)/log(10),'b-',tout2+st,quart2(3,:)/log(10),'r-',...
    mt0,lmtraj0,'b--',ft0,lftraj0,'r--',mt1,lmtraj1,'b-.',ft1,lftraj1,'r-.',...
    'LineWidth',1.5);
% set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
legend(h,{'Male, 3rd Quartile','Female, 3rd Quartile',...
    sprintf('Male Parameters, x0 = %2.1e',x0),...
    sprintf('Female Parameters, x0 = %2.1e',x0),...
    sprintf('Male Parameters, x0 = %2.1e',x1),...
    sprintf('Female Parameters, x0 = %2.1e',x1)},'Location','Best','FontSize',fontsizes(4));
xlabel('Days after Inoculation','FontSize',fontsizes(2));
ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
box off;
% title(sprintf('%s',desc));
% get detection limit;
% conv = getConverter();
% dlim = conv.a2v(3);
% dlim = 100;
if rnr == 1
    axlims = [st-0.5 st+tout1(end)+0.5 5.2 7.2];
elseif rnr == 0
    axlims = [st-0.5 st+tout1(end)+0.5 5.5 7.5];
else
    axlims = [st-0.5 st+tout1(end)+0.5 5.2 8];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesBySize(gdata,rnr,desc,tsel,mkey,fkey,size_range)
% get size trajectory of male and females
dt = 0.01;
if rnr == 0
    tf = 13;
    x0 = 3.0e5;
    x1 = 3.0e6;
else
    tf = 12;
    x0 = 7e5;
    x1 = 2e5;
end
[mtraj0,mt0] = getSizeTrajectory_v3_3(mkey,x0,dt,tf);
lmtraj0 = log10(mtraj0);
[ftraj0,ft0] = getSizeTrajectory_v3_3(fkey,x0,dt,tf);
lftraj0 = log10(ftraj0);
[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey,x1,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey,x1,dt,tf);
lftraj1 = log10(ftraj1);
% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .75 .85]);
fontsizes = [12 15 15 8]*2; 
if(rnr == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
elseif(rnr == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
end

% find tumors within size range
in_range1 = false(size(tlr1,1),1);%find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2))));
in_range2 = false(size(tlr2,1),1);%find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2))));

% let's scatter plot and quartile
split = [0.25 0.5 0.75]';
quart1 = zeros(length(split),length(tout1));
quart2 = zeros(length(split),length(tout2));
for i = 1:length(tout1)
%     mh1 = plot((tout1(i)+1-.05)*ones(1,avg_count1(i)),log_val1(1:avg_count1(i),i)/log(10),'bx');
%     hold on;
%     mh2 = plot((tout2(i)+1.05)*ones(1,avg_count2(i)),log_val2(1:avg_count2(i),i)/log(10),'r+');
    % get quartiles
    [sval1,~] = sort(log_val1(1:avg_count1(i),i));
    [sval2,~] = sort(log_val2(1:avg_count2(i),i));
    iq1_lo = 1+floor(split*(avg_count1(i)-1));
    iq1_hi = 1+ceil(split*(avg_count1(i)-1));
    iq2_lo = 1+floor(split*(avg_count2(i)-1));
    iq2_hi = 1+ceil(split*(avg_count2(i)-1));
    nq1_lo = (iq1_lo-1)/(avg_count1(i)-1);
    nq1_hi = (iq1_hi-1)/(avg_count1(i)-1);
    nq2_lo = (iq2_lo-1)/(avg_count2(i)-1);
    nq2_hi = (iq2_hi-1)/(avg_count2(i)-1);
    slp1 = (sval1(iq1_hi)-sval1(iq1_lo))./(nq1_hi-nq1_lo);
    slp1(nq1_hi-nq1_lo==0) = 0;
    quart1(:,i) = sval1(iq1_lo)+slp1.*(split-nq1_lo);
    slp2 = (sval2(iq2_hi)-sval2(iq2_lo))./(nq2_hi-nq2_lo);
    slp2(nq2_hi-nq2_lo==0) = 0;
    quart2(:,i) = sval2(iq2_lo)+slp2.*(split-nq2_lo);
end

h = plot(tout1+st,quart1(3,:)/log(10),'b-',tout2+st,quart2(3,:)/log(10),'r-',...
    mt0,lmtraj0,'b--',ft0,lftraj0,'r--',mt1,lmtraj1,'b-.',ft1,lftraj1,'r-.',...
    'LineWidth',1.5);
hold on;
for n = 1:length(in_range1)
    nn = in_range1(n);
    endi = find(tlr1(nn,:)~=0,1,'last');
    plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'bo-');
end
for n = 1:length(in_range2)
    nn = in_range2(n);
    endi = find(tlr2(nn,:)~=0,1,'last');
    plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'ro-');
end
hold off;
% set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
legend(h,{'Male, 3rd Quartile','Female, 3rd Quartile',...
    sprintf('Male Parameters, x0 = %2.1e',x0),...
    sprintf('Female Parameters, x0 = %2.1e',x0),...
    sprintf('Male Parameters, x0 = %2.1e',x1),...
    sprintf('Female Parameters, x0 = %2.1e',x1)},'Location','BestOutside',...
    'FontSize',fontsizes(4),'Box','off');
xlabel('Days after Inoculation','FontSize',fontsizes(2));
ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
box off;
set(gca,'FontSize',fontsizes(1));
% title(sprintf('%s',desc));
% get detection limit;
% conv = getConverter();
% dlim = conv.a2v(3);
% dlim = 100;
if rnr == 1
    axlims = [st-0.5 st+tout1(end)+0.5 5.2 7];
elseif rnr == 0
    axlims = [st-0.5 st+tout1(end)+0.5 5.5 7];
else
    axlims = [st-0.5 st+tout1(end)+0.5 5.2 8];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesBySize1(gdata,group_sel,title_text,tsel,mkey,fkey,size_range)
% plotting median and mean

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .75 .85]);
fontsizes = [6 8 8 7]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2))));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2))));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0m = getDay0Size(mkey,x1);
x0f = getDay0Size(fkey,x1);
        
% [mtraj0,mt0] = getSizeTrajectory_v3_3(mkey,x0,dt,tf);
% lmtraj0 = log10(mtraj0);
% [ftraj0,ft0] = getSizeTrajectory_v3_3(fkey,x0,dt,tf);
% lftraj0 = log10(ftraj0);
[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey,x0m,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey,x0f,dt,tf);
lftraj1 = log10(ftraj1);

h = plot(t1,med1,'b-',t2,med2,'r-',...
    t1,avg1,'b--',t2,avg2,'r--',...
    mt1,lmtraj1,'b-.',ft1,lftraj1,'r-.',...
    'LineWidth',1.5);
% unplotted
%     mt0,lmtraj0,'b.',ft0,lftraj0,'r.',...
%     sprintf('Male Parameters, x0 = %2.1e',x0),...
%     sprintf('Female Parameters, x0 = %2.1e',x0),...
    
% hold on;
% for n = 1:length(in_range1)
%     nn = in_range1(n);
%     endi = find(tlr1(nn,:)~=0,1,'last');
%     plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'b+-');
% end
% for n = 1:length(in_range2)
%     nn = in_range2(n);
%     endi = find(tlr2(nn,:)~=0,1,'last');
%     plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'r+-');
% end
% hold off;
% set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
legend(h,{sprintf('Male, Median, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range1)),...
    sprintf('Female, Median, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range2))...
    sprintf('Male, Average, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range1)),...
    sprintf('Female, Average, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range2))...
    sprintf('Male Parameters, x0 = %2.1e',x1),...
    sprintf('Female Parameters, x0 = %2.1e',x1)},'Location','Best','FontSize',fontsizes(4));
xlabel('Days after Inoculation','FontSize',fontsizes(2));
ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
title(title_text);
box off;
% title(sprintf('%s',desc));
% get detection limit;
% conv = getConverter();
% dlim = conv.a2v(3);
% dlim = 100;

max_val = ceil(max([med1(:);med2(:);lmtraj1(:);lftraj1(:)]));
if group_sel == 1
    axlims = [0 14 max_val-3 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-3 max_val];
else
    axlims = [0 14 max_val-3 max_val];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesBySize2(gdata,group_sel,title_text,tsel,mkey,fkey,size_range,gender)
% plotting median and mean along side individual trajectories

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .75 .85]);
fontsizes = [6 8 8 7]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2))));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2))));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0m = getDay0Size(mkey,x1);
x0f = getDay0Size(fkey,x1);
        
% [mtraj0,mt0] = getSizeTrajectory_v3_3(mkey,x0,dt,tf);
% lmtraj0 = log10(mtraj0);
% [ftraj0,ft0] = getSizeTrajectory_v3_3(fkey,x0,dt,tf);
% lftraj0 = log10(ftraj0);
[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey,x0m,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey,x0f,dt,tf);
lftraj1 = log10(ftraj1);

% different male and female plot
if(gender==0)
    h = plot(t1,med1,'b-',...
        t1,avg1,'b--',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'b+-');
    end
    hold off;
    % set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
    legend(h,{sprintf('Male, Median, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range1)),...
        sprintf('Male, Average, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range1))},...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(t2,med2,'r-',...
        t2,avg2,'r--',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'r+-');
    end
    hold off;
    % set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
    legend(h,{sprintf('Female, Median, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range2))...
        sprintf('Female, Average, %2.1e<x1<%2.1e, n=%i',size_range(1),size_range(2),length(in_range2))},...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4 max_val];
else
    axlims = [0 21 max_val-4 max_val];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesBySize3(gdata,group_sel,title_text,tsel,mkey,fkey,size_range,gender)
% plotting theory curve along side individual trajectories

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .75 .85]);
fontsizes = [6 8 8 7]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2))));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2))));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0m = getDay0Size(mkey,x1);
x0f = getDay0Size(fkey,x1);
        
% [mtraj0,mt0] = getSizeTrajectory_v3_3(mkey,x0,dt,tf);
% lmtraj0 = log10(mtraj0);
% [ftraj0,ft0] = getSizeTrajectory_v3_3(fkey,x0,dt,tf);
% lftraj0 = log10(ftraj0);
[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey,x0m,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey,x0f,dt,tf);
lftraj1 = log10(ftraj1);

% different male and female plot
if(gender==0)
    h = plot(mt1,lmtraj1,'b--',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'b+-');
    end
    hold off;
    % set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
    legend(h,{sprintf('Male Prediction, x0 = %2.1e',x1),'Male Tumor'},...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(ft1,lftraj1,'r--',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'r+-');
    end
    hold off;
    % set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
    legend(h,{sprintf('Female Prediction, x0 = %2.1e',x1),'Female Tumor'},...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4 max_val];
else
    axlims = [0 21 max_val-4 max_val];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesBySize4(gdata,group_sel,title_text,tsel,mkey1,fkey1,mkey2,fkey2,size_range,gender)
% plotting theory curve along side individual trajectories

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .65 .55]);
fontsizes = [6 8 8 7]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2))));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2))));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0m = getDay0Size(mkey1,x1);
x0f = getDay0Size(fkey1,x1);

[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey1,x0m,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey1,x0f,dt,tf);
lftraj1 = log10(ftraj1);

[mtraj2,mt2] = getSizeTrajectory_v3_3(mkey2,x0m,dt,tf);
lmtraj2 = log10(mtraj2);
[ftraj2,ft2] = getSizeTrajectory_v3_3(fkey2,x0f,dt,tf);
lftraj2 = log10(ftraj2);

% different male and female plot
if(gender==0)
    h = plot(mt1,lmtraj1,'k--',mt2,lmtraj2,'k.-',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'ko-');
    end
    hold off;
    % set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
    legend(h,{sprintf('Male Prediction 1, x0 = %2.1e',x1),...
        sprintf('Male Prediction 2, x0 = %2.1e',x1),'Male Tumor'},...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(ft1,lftraj1,'k--',ft2,lftraj2,'k.-',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'ko-');
    end
    hold off;
    % set(h,{'LineStyle'},{'-.';'-';'--';'-.';'-';'--'});
    legend(h,{sprintf('Female Prediction 1, x0 = %2.1e',x1),...
        sprintf('Female Prediction 2, x0 = %2.1e',x1),'Female Tumor'},...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4 max_val];
else
    axlims = [0 21 max_val-4 max_val];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesBySize5(gdata,group_sel,title_text,tsel,mkey1,fkey1,mkey2,fkey2,size_range,gender)
% plotting theory curve along side individual trajectories
% sort tumor curves by ages and color differently

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('position',[769.8000 41.8000 766.4000 740.8000]);
fontsizes = [9 12 12 10]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
    age1 = gdata(1).age_rad(gdata(1).viable_rad);
    age2 = gdata(2).age_rad(gdata(2).viable_rad);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
    age1 = gdata(1).age_norad(gdata(1).viable_norad);
    age2 = gdata(2).age_norad(gdata(2).viable_norad);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
    age1 = gdata(1).age_norad2(gdata(1).viable_norad2);
    age2 = gdata(2).age_norad2(gdata(2).viable_norad2);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2)))&(age1>300));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2)))&(age2>300));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0m = getDay0Size(mkey1,x1);
x0f = getDay0Size(fkey1,x1);

[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey1,x0m,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey1,x0f,dt,tf);
lftraj1 = log10(ftraj1);

[mtraj2,mt2] = getSizeTrajectory_v3_3(mkey2,x0m,dt,tf);
lmtraj2 = log10(mtraj2);
[ftraj2,ft2] = getSizeTrajectory_v3_3(fkey2,x0f,dt,tf);
lftraj2 = log10(ftraj2);

% change colors based on young and old fish
age_cutoff = 300;
age_crit1 = (age1>age_cutoff)+1;
age_crit2 = (age2>age_cutoff)+1;
age_colors = [0 0.447 0.741; 0.850 0.325 0.098];
age_markers = {'s','o'};

% different male and female plot
if(gender==0)
    h = plot(mt1,lmtraj1,'k--',mt2,lmtraj2,'k-.',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        hn = plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit1(nn)},...
            'Color',age_colors(age_crit1(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,2+length(in_range1));
    leg_text(1:2) = {sprintf('Male Prediction 1, x0 = %2.1e',x1),...
        sprintf('Male Prediction 2, x0 = %2.1e',x1)};
    first_young = find(ismember(in_range1,find(age_crit1==1)),1,'first');
    first_old = find(ismember(in_range1,find(age_crit1==2)),1,'first');
    for ii = 3:length(leg_text)
        if(ii==first_young+2)
            leg_text{first_young+2} = 'Male, 6 months old';
        elseif(ii==first_old+2)
            leg_text{first_old+2} = 'Male, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1:2 first_young+2 first_old+2]),...
        leg_text([1:2 first_young+2 first_old+2]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(ft1,lftraj1,'k--',ft2,lftraj2,'k-.',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        hn = plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit2(nn)},...
            'Color',age_colors(age_crit2(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,2+length(in_range2));
    leg_text(1:2) = {sprintf('Female Prediction 1, x0 = %2.1e',x1),...
        sprintf('Female Prediction 2, x0 = %2.1e',x1)};
    first_young = find(ismember(in_range2,find(age_crit2==1)),1,'first');
    first_old = find(ismember(in_range2,find(age_crit2==2)),1,'first');
    for ii = 3:length(leg_text)
        if(ii==first_young+2)
            leg_text{first_young+2} = 'Female, 6 months old';
        elseif(ii==first_old+2)
            leg_text{first_old+2} = 'Female, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1:2 first_young+2 first_old+2]),...
        leg_text([1:2 first_young+2 first_old+2]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4.0 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4.0 max_val];
else
    axlims = [0 21 max_val-4.0 max_val];
end
axis(axlims);

legend boxoff;

set(gca,'FontSize',fontsizes(1));

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesAreaError(gdata,group_sel,title_text,tsel,mid_key,hi_key,lo_key,size_range,gender)
% plotting theory curve along side individual trajectories
% sort tumor curves by ages and color differently

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('position',[769.8000 41.8000 766.4000 740.8000]);
fontsizes = [9 12 12 10]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
    age1 = gdata(1).age_rad(gdata(1).viable_rad);
    age2 = gdata(2).age_rad(gdata(2).viable_rad);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
    age1 = gdata(1).age_norad(gdata(1).viable_norad);
    age2 = gdata(2).age_norad(gdata(2).viable_norad);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
    age1 = gdata(1).age_norad2(gdata(1).viable_norad2);
    age2 = gdata(2).age_norad2(gdata(2).viable_norad2);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2)))&(age1>300));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2)))&(age2>300));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0 = getDay0Size(mid_key,x1);

[mid_traj,mid_t] = getSizeTrajectory_v3_3(mid_key,x0,dt,tf);
lmid_traj = log10(mid_traj);

[hi_traj,hi_t] = getSizeTrajectory_v3_3(hi_key,x0,dt,tf);
lhi_traj = log10(hi_traj);

[lo_traj,lo_t] = getSizeTrajectory_v3_3(lo_key,x0,dt,tf);
llo_traj = log10(lo_traj);

% change colors based on young and old fish
age_cutoff = 300;
age_crit1 = (age1>age_cutoff)+1;
age_crit2 = (age2>age_cutoff)+1;
age_colors = [0 0.447 0.741; 0.850 0.325 0.098];
age_markers = {'s','o'};

% different male and female plot
if(gender==0)
    h = plot(mid_t,lmid_traj,'k--','LineWidth',1.5);
    hold on;
    
    fill_x = [mid_t fliplr(mid_t)];
    fill_y = [lhi_traj fliplr(max(llo_traj,0))];
    fill(fill_x,fill_y,'g','faceAlpha',0.33,'edgecolor','none');
    
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        hn = plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit1(nn)},...
            'Color',age_colors(age_crit1(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,1+length(in_range1));
    leg_text{1} = sprintf('Male Prediction, x0 = %2.1e',x1);
    first_young = find(ismember(in_range1,find(age_crit1==1)),1,'first');
    first_old = find(ismember(in_range1,find(age_crit1==2)),1,'first');
    for ii = 2:length(leg_text)
        if(ii==first_young+1)
            leg_text{first_young+1} = 'Male, 6 months old';
        elseif(ii==first_old+1)
            leg_text{first_old+1} = 'Male, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1 first_young+1 first_old+1]),...
        leg_text([1 first_young+1 first_old+1]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(mid_t,lmid_traj,'k--','LineWidth',1.5);
    hold on;
    
    fill_x = [mid_t fliplr(mid_t)];
    fill_y = [lhi_traj fliplr(max(llo_traj,0))];
    fill(fill_x,fill_y,'g','faceAlpha',0.33,'edgecolor','none');
    
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        hn = plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit2(nn)},...
            'Color',age_colors(age_crit2(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,1+length(in_range2));
    leg_text{1} = sprintf('Female Prediction, x0 = %2.1e',x1);
    first_young = find(ismember(in_range2,find(age_crit2==1)),1,'first');
    first_old = find(ismember(in_range2,find(age_crit2==2)),1,'first');
    for ii = 2:length(leg_text)
        if(ii==first_young+1)
            leg_text{first_young+1} = 'Female, 6 months old';
        elseif(ii==first_old+1)
            leg_text{first_old+1} = 'Female, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1 first_young+1 first_old+1]),...
        leg_text([1 first_young+1 first_old+1]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4.0 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4.0 max_val];
else
    axlims = [0 21 max_val-4.0 max_val];
end
axis(axlims);

legend boxoff;

set(gca,'FontSize',fontsizes(1));

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesAreaError2(gdata,group_sel,title_text,tsel,hilo,size_range,gender)
% plotting theory curve along side individual trajectories
% sort tumor curves by ages and color differently

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('position',[769.8000 41.8000 766.4000 740.8000]);
fontsizes = [9 12 12 10]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
    age1 = gdata(1).age_rad(gdata(1).viable_rad);
    age2 = gdata(2).age_rad(gdata(2).viable_rad);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
    age1 = gdata(1).age_norad(gdata(1).viable_norad);
    age2 = gdata(2).age_norad(gdata(2).viable_norad);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
    age1 = gdata(1).age_norad2(gdata(1).viable_norad2);
    age2 = gdata(2).age_norad2(gdata(2).viable_norad2);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2)))&(age1>300));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2)))&(age2>300));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
mid_key = hilo{1}.key;
x0 = getDay0Size(mid_key,x1);

[mid_traj,mid_t] = getSizeTrajectory_v3_3(mid_key,x0,dt,tf);
lmid_traj = log10(mid_traj);

hilo_trajs = zeros(length(hilo)-1,length(mid_traj));

for hh = 1:size(hilo_trajs,1)
    [hilo_trajs(hh,:),~] = getSizeTrajectory_v3_3(hilo{hh+1}.key,x0,dt,tf);
end

hi_traj = max(hilo_trajs,[],1);
lhi_traj = log10(hi_traj);

lo_traj = min(hilo_trajs,[],1);
llo_traj = log10(lo_traj);

% change colors based on young and old fish
age_cutoff = 300;
age_crit1 = (age1>age_cutoff)+1;
age_crit2 = (age2>age_cutoff)+1;
age_colors = [0 0.447 0.741; 0.850 0.325 0.098];
age_markers = {'s','o'};

% different male and female plot
if(gender==0)
    h = plot(mid_t,lmid_traj,'k--','LineWidth',1.5);
    hold on;
    
    fill_x = [mid_t fliplr(mid_t)];
    fill_y = [lhi_traj fliplr(max(llo_traj,0))];
    fill(fill_x,fill_y,'g','faceAlpha',0.33,'edgecolor','none');
    
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        hn = plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit1(nn)},...
            'Color',age_colors(age_crit1(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,1+length(in_range1));
    leg_text{1} = sprintf('Male Prediction, x0 = %2.1e',x1);
    first_young = find(ismember(in_range1,find(age_crit1==1)),1,'first');
    first_old = find(ismember(in_range1,find(age_crit1==2)),1,'first');
    for ii = 2:length(leg_text)
        if(ii==first_young+1)
            leg_text{first_young+1} = 'Male, 6 months old';
        elseif(ii==first_old+1)
            leg_text{first_old+1} = 'Male, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1 first_young+1 first_old+1]),...
        leg_text([1 first_young+1 first_old+1]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(mid_t,lmid_traj,'k--','LineWidth',1.5);
    hold on;
    
    fill_x = [mid_t fliplr(mid_t)];
    fill_y = [lhi_traj fliplr(max(llo_traj,0))];
    fill(fill_x,fill_y,'g','faceAlpha',0.33,'edgecolor','none');
    
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        hn = plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit2(nn)},...
            'Color',age_colors(age_crit2(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,1+length(in_range2));
    leg_text{1} = sprintf('Female Prediction, x0 = %2.1e',x1);
    first_young = find(ismember(in_range2,find(age_crit2==1)),1,'first');
    first_old = find(ismember(in_range2,find(age_crit2==2)),1,'first');
    for ii = 2:length(leg_text)
        if(ii==first_young+1)
            leg_text{first_young+1} = 'Female, 6 months old';
        elseif(ii==first_old+1)
            leg_text{first_old+1} = 'Female, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1 first_young+1 first_old+1]),...
        leg_text([1 first_young+1 first_old+1]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4.0 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4.0 max_val];
else
    axlims = [0 21 max_val-4.0 max_val];
end
axis(axlims);

legend boxoff;

set(gca,'FontSize',fontsizes(1));

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotWithTrajectoriesAreaError3(gdata,group_sel,title_text,tsel,hilo1,hilo2,size_range,gender)
% plotting theory curve along side individual trajectories
% sort tumor curves by ages and color differently

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('position',[769.8000 41.8000 766.4000 740.8000]);
fontsizes = [9 12 12 10]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
    age1 = gdata(1).age_rad(gdata(1).viable_rad);
    age2 = gdata(2).age_rad(gdata(2).viable_rad);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
    age1 = gdata(1).age_norad(gdata(1).viable_norad);
    age2 = gdata(2).age_norad(gdata(2).viable_norad);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
    age1 = gdata(1).age_norad2(gdata(1).viable_norad2);
    age2 = gdata(2).age_norad2(gdata(2).viable_norad2);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2)))&(age1>300));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2)))&(age2>300));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
mid_key = hilo1{1}.key;
x0 = getDay0Size(mid_key,x1);

[mid_traj1,mid_t1] = getSizeTrajectory_v3_3(mid_key,x0,dt,tf);
lmid_traj1 = log10(mid_traj1);

hilo_trajs = zeros(length(hilo1)-1,length(mid_traj1));

for hh = 1:size(hilo_trajs,1)
    [hilo_trajs(hh,:),~] = getSizeTrajectory_v3_3(hilo1{hh+1}.key,x0,dt,tf);
end

hi_traj = max(hilo_trajs,[],1);
lhi_traj1 = log10(hi_traj);

lo_traj = min(hilo_trajs,[],1);
llo_traj1 = log10(lo_traj);

% second time with second set of parameters
mid_key = hilo2{1}.key;
x0 = getDay0Size(mid_key,x1);

[mid_traj2,mid_t2] = getSizeTrajectory_v3_3(mid_key,x0,dt,tf);
lmid_traj2 = log10(mid_traj2);

hilo_trajs = zeros(length(hilo2)-1,length(mid_traj1));

for hh = 1:size(hilo_trajs,1)
    [hilo_trajs(hh,:),~] = getSizeTrajectory_v3_3(hilo2{hh+1}.key,x0,dt,tf);
end

hi_traj = max(hilo_trajs,[],1);
lhi_traj2 = log10(hi_traj);

lo_traj = min(hilo_trajs,[],1);
llo_traj2 = log10(lo_traj);


% change colors based on young and old fish
age_cutoff = 300;
age_crit1 = (age1>age_cutoff)+1;
age_crit2 = (age2>age_cutoff)+1;
age_colors = [0 0.447 0.741; 0.850 0.325 0.098];
age_markers = {'s','o'};

% different male and female plot
if(gender==0)
    h = plot(mid_t1,lmid_traj1,'g--',mid_t2,lmid_traj2,'b--','LineWidth',1.5);
    hold on;
    
    fill_x1 = [mid_t1 fliplr(mid_t1)];
    fill_y1 = max([lhi_traj1 fliplr(llo_traj1)],0);
    fill(fill_x1,fill_y1,'g','faceAlpha',0.33,'edgecolor','none');
    
    fill_x2 = [mid_t2 fliplr(mid_t2)];
    fill_y2 = max([lhi_traj2 fliplr(llo_traj2)],0);
    fill(fill_x2,fill_y2,'b','faceAlpha',0.33,'edgecolor','none');
    
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        hn = plot(tlr1(nn,1:endi),log_rad1(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit1(nn)},...
            'Color',age_colors(age_crit1(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,2+length(in_range1));
    leg_text{1} = sprintf('Prediction, Day 4 immunity, x0 = %2.1e',x1);
    leg_text{2} = sprintf('Prediction, Day 11 immunity, x0 = %2.1e',x1);
    first_young = find(ismember(in_range1,find(age_crit1==1)),1,'first');
    first_old = find(ismember(in_range1,find(age_crit1==2)),1,'first');
    for ii = 3:length(leg_text)
        if(ii==first_young+2)
            leg_text{first_young+2} = 'Male, 6 months old';
        elseif(ii==first_old+2)
            leg_text{first_old+2} = 'Male, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1:2 first_young+2 first_old+2]),...
        leg_text([1:2 first_young+2 first_old+2]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max([med1(:); log_rad1(:)/log(10)]));
else
    h = plot(mid_t1,lmid_traj1,'g--',mid_t2,lmid_traj2,'b--','LineWidth',1.5);
    hold on;
    
    fill_x1 = [mid_t1 fliplr(mid_t1)];
    fill_y1 = max([lhi_traj1 fliplr(llo_traj1)],0);
    fill(fill_x1,fill_y1,'g','faceAlpha',0.33,'edgecolor','none');
    
    fill_x2 = [mid_t2 fliplr(mid_t2)];
    fill_y2 = max([lhi_traj2 fliplr(llo_traj2)],0);
    fill(fill_x2,fill_y2,'b','faceAlpha',0.33,'edgecolor','none');
    
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        hn = plot(tlr2(nn,1:endi),log_rad2(nn,1:endi)/log(10),'-',...
            'Marker',age_markers{age_crit2(nn)},...
            'Color',age_colors(age_crit2(nn),:),...
            'LineWidth',1.25,'MarkerSize',7);
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,2+length(in_range2));
    leg_text{1} = sprintf('Prediction, Day 4 immunity, x0 = %2.1e',x1);
    leg_text{2} = sprintf('Prediction, Day 11 immunity, x0 = %2.1e',x1);
    first_young = find(ismember(in_range2,find(age_crit2==1)),1,'first');
    first_old = find(ismember(in_range2,find(age_crit2==2)),1,'first');
    for ii = 3:length(leg_text)
        if(ii==first_young+2)
            leg_text{first_young+2} = 'Female, 6 months old';
        elseif(ii==first_old+2)
            leg_text{first_old+2} = 'Female, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1:2 first_young+2 first_old+2]),...
        leg_text([1:2 first_young+2 first_old+2]),...
        'Location','SouthWest','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max([med2(:); log_rad2(:)/log(10)]));
end

if group_sel == 1
    axlims = [0 14 max_val-4.0 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-4.0 max_val];
else
    axlims = [0 21 max_val-4.0 max_val];
end
axis(axlims);

legend boxoff;

set(gca,'FontSize',fontsizes(1));

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [] = plotRelativeSize(gdata,group_sel,title_text,tsel,mkey1,fkey1,mkey2,fkey2,size_range,gender)
% plotting theory curve along side individual trajectories
% sort tumor curves by ages and color differently

% plot male and female together
% rnr = 0 for rad and 1 for norad 
figure('units','normalized','outerposition',[.05 .05 .65 .55]);
fontsizes = [6 8 8 7]*2; 
if(group_sel == 0)
    st = gdata(1).trad(1,1);
    tout1 = gdata(1).tout_r(tsel);
    log_val1 = gdata(1).log_val_r(:,tsel);
    avg_count1 = gdata(1).avg_count_r(tsel);
    log_avg1 = gdata(1).log_avg_r(tsel);
    log_med1 = gdata(1).log_med_r(tsel);
    tout2 = gdata(2).tout_r(tsel);
    log_val2 = gdata(2).log_val_r(:,tsel);
    avg_count2 = gdata(2).avg_count_r(tsel);
    log_avg2 = gdata(2).log_avg_r(tsel);
    log_med2 = gdata(2).log_med_r(tsel);
    tlr1 = gdata(1).trad(:,tsel);
    tlr2 = gdata(2).trad(:,tsel);
    log_rad1 = gdata(1).log_rad(:,tsel);
    log_rad2 = gdata(2).log_rad(:,tsel);
    age1 = gdata(1).age_rad(gdata(1).viable_rad);
    age2 = gdata(2).age_rad(gdata(2).viable_rad);
elseif(group_sel == 1)
    st = gdata(1).tnorad(1,1);
    tout1 = gdata(1).tout_nr(tsel);
    log_val1 = gdata(1).log_val_nr(:,tsel);
    avg_count1 = gdata(1).avg_count_nr(tsel);
    log_avg1 = gdata(1).log_avg_nr(tsel);
    log_med1 = gdata(1).log_med_nr(tsel);
    tout2 = gdata(2).tout_nr(tsel);
    log_val2 = gdata(2).log_val_nr(:,tsel);
    avg_count2 = gdata(2).avg_count_nr(tsel);
    log_avg2 = gdata(2).log_avg_nr(tsel);
    log_med2 = gdata(2).log_med_nr(tsel);
    tlr1 = gdata(1).tnorad(:,tsel);
    tlr2 = gdata(2).tnorad(:,tsel);
    log_rad1 = gdata(1).log_norad(:,tsel);
    log_rad2 = gdata(2).log_norad(:,tsel);
    age1 = gdata(1).age_norad(gdata(1).viable_norad);
    age2 = gdata(2).age_norad(gdata(2).viable_norad);
else
    st = gdata(1).tnorad2(1,1);
    tout1 = gdata(1).tout_nr2(tsel);
    log_val1 = gdata(1).log_val_nr2(:,tsel);
    avg_count1 = gdata(1).avg_count_nr2(tsel);
    log_avg1 = gdata(1).log_avg_nr2(tsel);
    log_med1 = gdata(1).log_med_nr2(tsel);
    tout2 = gdata(2).tout_nr2(tsel);
    log_val2 = gdata(2).log_val_nr2(:,tsel);
    avg_count2 = gdata(2).avg_count_nr2(tsel);
    log_avg2 = gdata(2).log_avg_nr2(tsel);
    log_med2 = gdata(2).log_med_nr2(tsel);
    tlr1 = gdata(1).tnorad2(:,tsel);
    tlr2 = gdata(2).tnorad2(:,tsel);
    log_rad1 = gdata(1).log_norad2(:,tsel);
    log_rad2 = gdata(2).log_norad2(:,tsel);
    age1 = gdata(1).age_norad2(gdata(1).viable_norad2);
    age2 = gdata(2).age_norad2(gdata(2).viable_norad2);
end

% find tumors within size range
in_range1 = find((tlr1(:,1)==1)&(log_rad1(:,1)>=log(size_range(1)))&(log_rad1(:,1)<=log(size_range(2)))&(age1>300));
in_range2 = find((tlr2(:,1)==1)&(log_rad2(:,1)>=log(size_range(1)))&(log_rad2(:,1)<=log(size_range(2)))&(age2>300));

% get medians and averages for values in range (THIS ASSUMES TIME POINTS
% MATCH WHICH MAY NOT ALWAYS BE TRUE)
if(~isempty(in_range1))
    t1 = tlr1(in_range1(1),:);
    med1 = median(log_rad1(in_range1,:)/log(10),1);
    avg1 = mean(log_rad1(in_range1,:)/log(10),1);
else
    t1 = [1];
    med1 = log10(2e5);
    avg1 = log10(2e5);
end
if(~isempty(in_range2))
    t2 = tlr2(in_range2(1),:);
    med2 = median(log_rad2(in_range2,:)/log(10),1);
    avg2 = mean(log_rad2(in_range2,:)/log(10),1);
else
    t2 = [1];
    med2 = log10(2e5);
    avg2 = log10(2e5);
end

% get size trajectory of male and females
dt = 0.01;
if group_sel == 0
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
else
    tf = 21;
%     x0 = size_range(1);
    x1 = 10^(0.5*(med1(1)+med2(1)));
end

% get proper starting sizes to match Day 1 size
x0m = getDay0Size(mkey1,x1);
x0f = getDay0Size(fkey1,x1);

[mtraj1,mt1] = getSizeTrajectory_v3_3(mkey1,x0m,dt,tf);
lmtraj1 = log10(mtraj1);
[ftraj1,ft1] = getSizeTrajectory_v3_3(fkey1,x0f,dt,tf);
lftraj1 = log10(ftraj1);

[mtraj2,mt2] = getSizeTrajectory_v3_3(mkey2,x0m,dt,tf);
lmtraj2 = log10(mtraj2);
[ftraj2,ft2] = getSizeTrajectory_v3_3(fkey2,x0f,dt,tf);
lftraj2 = log10(ftraj2);

% change colors based on young and old fish
age_cutoff = 300;
age_crit1 = (age1>age_cutoff)+1;
age_crit2 = (age2>age_cutoff)+1;
age_colors = [0 0.447 0.741; 0.850 0.325 0.098];
age_markers = {'s','o'};

% different male and female plot
if(gender==0)
    h = plot(mt1,lmtraj1-log10(x1),'k--',mt2,lmtraj2-log10(x1),'k.-',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range1)
        nn = in_range1(n);
        endi = find(tlr1(nn,:)~=0,1,'last');
        hn = plot(tlr1(nn,1:endi),(log_rad1(nn,1:endi)-log_rad1(nn,1))/log(10),'-',...
            'Marker',age_markers{age_crit1(nn)},...
            'Color',age_colors(age_crit1(nn),:));
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,2+length(in_range1));
    leg_text(1:2) = {sprintf('Male Prediction 1, x0 = %2.1e',x1),...
        sprintf('Male Prediction 2, x0 = %2.1e',x1)};
    first_young = find(ismember(in_range1,find(age_crit1==1)),1,'first');
    first_old = find(ismember(in_range1,find(age_crit1==2)),1,'first');
    for ii = 3:length(leg_text)
        if(ii==first_young+2)
            leg_text{first_young+2} = 'Male, 6 months old';
        elseif(ii==first_old+2)
            leg_text{first_old+2} = 'Male, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1:2 first_young+2 first_old+2]),...
        leg_text([1:2 first_young+2 first_old+2]),...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    max_val = ceil(max(max((log_rad1-log_rad1(:,1))/log(10))));
else
    h = plot(ft1,lftraj1-log10(x1),'k--',ft2,lftraj2-log10(x1),'k.-',...
        'LineWidth',1.5);
    
    hold on;
    for n = 1:length(in_range2)
        nn = in_range2(n);
        endi = find(tlr2(nn,:)~=0,1,'last');
        hn = plot(tlr2(nn,1:endi),(log_rad2(nn,1:endi)-log_rad2(nn,1))/log(10),'-',...
            'Marker',age_markers{age_crit2(nn)},...
            'Color',age_colors(age_crit2(nn),:));
        h = [h;hn];
    end
    hold off;
    
    leg_text = cell(1,2+length(in_range2));
    leg_text(1:2) = {sprintf('Female Prediction 1, x0 = %2.1e',x1),...
        sprintf('Female Prediction 2, x0 = %2.1e',x1)};
    first_young = find(ismember(in_range2,find(age_crit2==1)),1,'first');
    first_old = find(ismember(in_range2,find(age_crit2==2)),1,'first');
    for ii = 3:length(leg_text)
        if(ii==first_young+2)
            leg_text{first_young+2} = 'Female, 6 months old';
        elseif(ii==first_old+2)
            leg_text{first_old+2} = 'Female, 20 months old';
        else
            leg_text{ii} = '';
        end
    end
    legend(h([1:2 first_young+2 first_old+2]),...
        leg_text([1:2 first_young+2 first_old+2]),...
        'Location','Best','FontSize',fontsizes(4));
    xlabel('Days after Inoculation','FontSize',fontsizes(2));
    ylabel('Log_{10} Tumor Size','FontSize',fontsizes(3));
    title(title_text);
    box off;
    
    max_val = ceil(max(max((log_rad2-log_rad2(:,1))/log(10))));
end

if group_sel == 1
    axlims = [0 14 max_val-8 max_val];
elseif group_sel == 0
    axlims = [0 14 max_val-8 max_val];
else
    axlims = [0 21 max_val-8 max_val];
end
axis(axlims);

% edges = linspace(axlims(3),axlims(4),51);
% ax1 = gca;
% axslp = ax1.Position(3)/(axlims(2)-axlims(1));
% for i = 1:length(tout1)
%     pwid = 0.025;
%     ax2pos = [ax1.Position(1)+axslp*(.5+tout1(i))-pwid ax1.Position(2) pwid ax1.Position(4)];
%     ax2 = axes('Position',ax2pos);
%     sval1 = histcounts(log_val1(1:avg_count1(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),-sval1,'b','EdgeColor','none','BarWidth',0.8);
%     axis([-10 0 axlims(3:4)]);
%     ax2.XTickLabel = [];
%     ax2.XTick = [];
%     ax2.YTickLabel = [];
%     ax2.YTick = [];
%     ax2.Visible = 'off';
%     ax2.XColor = 'none';
%     ax3pos = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) ax2.Position(3) ax2.Position(4)];
%     ax3 = axes('Position',ax3pos);
%     sval2 = histcounts(log_val2(1:avg_count2(i),i)/log(10),edges);
%     barh(0.5*(edges(2:end)+edges(1:end-1)),sval2,'r','EdgeColor','none','BarWidth',0.8);
%     axis([0 10 axlims(3:4)]);
%     ax3.XTickLabel = [];
%     ax3.XTick = [];
%     ax3.YTickLabel = [];
%     ax3.YTick = [];
%     ax3.Visible = 'off';
%     ax3.XColor = 'none';
% end
end

function [inj_locs] = assignLoc(locs,lab,desc)
% use tumor identifying info in lab and desc and location info in locs to mark inj
% location of each tumor

inj_locs = cell(size(lab,1),1);
% find correct tumor_id for each tumor in lab
new_fish_pos = lab(:,2)==1; % number of new fish = size(desc,1)
id = 0;
for i = 1:length(inj_locs)
    if(new_fish_pos(i))
        id = id+1;
    end
    inj_locs{i} = locs{desc(id,1)};
end

end

function [ages] = assignAge(fish_ages,lab,desc)
% use tumor identifying info in lab and desc and location info in locs to mark inj
% location of each tumor

ages = zeros(size(lab,1),1);
% find correct tumor_id for each tumor in lab
new_fish_pos = lab(:,2)==1; % number of new fish = size(desc,1)
id = 0;
for i = 1:length(ages)
    if(new_fish_pos(i))
        id = id+1;
    end
    ages(i) = fish_ages(desc(id,1));
end

end

function [x0] = getDay0Size(key,x1)
    a = key.GROWTH_PARAMETER1;
    b = key.GROWTH_PARAMETER2;
    x0 = (x1.^(1-b)-a*(1-b)*1).^(1/(1-b));
end