%% runCreateFishStruct_v2_4
%  Version 2.4
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 2/15/20
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Version History
%  1.4: works with CreateFishStruct_v3_1 (step2_v1_3)
%  1.5: steps to finalize results
%  2.1: poincare version
%  2.2: clustering_threshold_day1 set to 10 (calculated only to 25 days)

% clear variables;
baseName = 'ct10_lt015_ht04_021620_nosat';
fish_dir = 'C:\Users\David\Documents\Yinka\fish_images\';
max_tp = 20;
complete_times = cell(1,max_tp);
% complete_times = {'1','3','5','7','9','11','13','15',...}; 
for ti = 1:max_tp
    complete_times{ti} = num2str(2*ti-1);
end
% complete_times = {'0','1'};
% timePoint = {'1','3','5','7','9'};%,'11'}; 
% exposure times ([BF;GFP;RFP;COLOR])
B0130_exp_times = [1.6 1.6 1.6 1.6 1.6 1.6 1.6 1.6 1.6 1.6 1.6 1.6 1.6;...
                   4000 3000 3000 206.4 300 300 300 1000 666.4 666.4 666.4 300 1000;...
                   4000 6000 6000 6000 3000 3000 3000 3000 3000 3000 3000 3000 3000;...
                   2.8 2.8 2.8 2.8 2.8 2.8 2.8 2.8 2.8 2.8 2.8 4.2 2.8];
B0130R_max_tp = 13;
B0130NR_max_tp = 6;
B0213_exp_times = [1.6 1.6 1.6 1.6 1.6 1.6;...
                   1500 1000 666.4 666.4 300 1000;...
                   3000 3000 3000 3000 3000 3000;...
                   2.8 2.8 2.8 2.8 4.2 2.8];
B0213R_max_tp = 6;
B0213NR_max_tp = 6;
B0709_exp_times = [3 3 3 3 3 3 3 3 3;...
                   1000 1000 300 300 300 300 300 300 300;...
                   1000 1000 700 700 700 700 700 700 700;...
                   4 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B0709R_max_tp = 9;
B0709R2_max_tp = 4;
B0728_exp_times = [3 3 3 3 3;...
                   700 400 300 300 700;...
                   700 700 700 700 700;...
                   4.2 4.2 4.2 4.2 4.2];
B0728NR_max_tp = 5;
B0814_exp_times = [3 3 3 3 3 3 3 3;...
                   300 300 300 300 300 300 300 300;...
                   700 700 700 700 700 700 700 700;...
                   4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B0814R_max_tp = 8;
B0903_exp_times = [4 4 4 4 4 4 4 4;...
                   400 400 400 400 400 400 400 400;...
                   1000 1000 1000 1000 1000 1000 1000 1000;...
                   4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B0903NR_max_tp = 8;

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

B090619_exp_times = [4 4 4 4 4 3 3
                     700 700 700 550 550 550 550
                     700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B090619_max_tp = length(B090619_times);

B091119_exp_times = [4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
                     700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700
                     700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B091119_max_tp = length(B091119_times);

B091919_exp_times = [3 3 3 3 3 3 3 3 3 3 3
                     550 550 550 300 300 300 300 500 500 550 700
                     700 700 700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B091919_max_tp = length(B091919_times);

B092319_exp_times = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
                     550 400 250 250 175 175 175 250 300 300 350 400 350 400 300 400 350 350 350 350 400 400 400
                     700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2]; 
B092319_max_tp = length(B092319_times);

B092719_exp_times = [3 3 3 3 3 3
                     500 250 250 300 700 700
                     700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2];
B092719_max_tp = length(B092719_times);

B100119_exp_times = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
                     550 550 350 350 350 300 350 400 500 400 400 300 300 300 400 400 500 500 700
                     700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B100119_max_tp = length(B100119_times);

B100919_exp_times = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
                     600 350 350 400 500 550 700 700 700 700 700 700 700 700 700
                     700 700 700 700 700 700 700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B100919_max_tp = length(B100919_times);

B101319_exp_times = [3 3 3 3 3 3 3 3 3
                     300 250 300 200 200 200 700 700 700
                     700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B101319_max_tp = length(B101319_times);

B101719_exp_times = [3 3 3 3 3 3 3 3 3 3 3
                     500 400 250 200 175 200 200 200 300 300 300
                     700 700 700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B101719_max_tp = length(B101719_times);

B102119_exp_times = [3 3 3 3 3 3 3 3 3
                     300 300 250 200 200 300 400 400 700
                     700 700 700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B102119_max_tp = length(B102119_times);

B102519_exp_times = [3 3 3 3 3 3 3
                     400 300 250 300 250 300 400
                     700 700 700 700 700 700 700
                     4.2 4.2 4.2 4.2 4.2 4.2 4.2];
B102519_max_tp = length(B102519_times);

normal_exposure_times = [2; 300; 300; 3];

% picture naming convention
B0130_img_naming_conv = struct('color','.','brightfield','_TL Brightfield_3.',...
    'gfp','_Alexa Fluor 488_1.','rfp','_Rhodamine_2.','overlay','_c1-3.')';
B0709_img_naming_conv = struct('color','.','brightfield','_TL Brightfield_1.',...
    'gfp','_EGFP_2.','rfp','_mCherry_3.','overlay','_c1-3.')';

batch_count = 0;

B0130_RAD_baseDir = cell(1,B0130R_max_tp);
for ti = 1:B0130R_max_tp
    B0130_RAD_baseDir{ti} = [fish_dir 'td013017 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD'];
end
batch_count = batch_count+1;

B0130_NORAD_baseDir = cell(1,B0130NR_max_tp);
for ti = 1:B0130NR_max_tp
    B0130_NORAD_baseDir{ti} = [fish_dir 'td013017 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B0213_RAD_baseDir = cell(1,B0213R_max_tp);
for ti = 1:B0213R_max_tp
    B0213_RAD_baseDir{ti} = [fish_dir 'td021317 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD'];
end
batch_count = batch_count+1;

B0213_NORAD_baseDir = cell(1,B0213NR_max_tp);
for ti = 1:B0213NR_max_tp
    B0213_NORAD_baseDir{ti} = [fish_dir 'td021317 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B0709_RAD_baseDir = cell(1,B0709R_max_tp);
for ti = 1:B0709R_max_tp
    B0709_RAD_baseDir{ti} = [fish_dir 'td070917 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD'];
end
batch_count = batch_count+1;

B0709_RAD2_baseDir = cell(1,B0709R2_max_tp);
for ti = 1:B0709R2_max_tp
    B0709_RAD2_baseDir{ti} = [fish_dir 'td070917 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD2'];
end
batch_count = batch_count+1;

B0728_NORAD1_baseDir = cell(1,B0728NR_max_tp);
for ti = 1:B0728NR_max_tp
    B0728_NORAD1_baseDir{ti} = [fish_dir 'td072817 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/NO RAD1'];
end
batch_count = batch_count+1;

B0728_NORAD2_baseDir = cell(1,B0728NR_max_tp);
for ti = 1:B0728NR_max_tp
    B0728_NORAD2_baseDir{ti} = [fish_dir 'td072817 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/NO RAD2'];
end
batch_count = batch_count+1;

B0728_NORAD3_baseDir = cell(1,B0728NR_max_tp);
for ti = 1:B0728NR_max_tp
    B0728_NORAD3_baseDir{ti} = [fish_dir 'td072817 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/NO RAD3'];
end
batch_count = batch_count+1;

B0728_NORAD4_baseDir = cell(1,B0728NR_max_tp);
for ti = 1:B0728NR_max_tp
    B0728_NORAD4_baseDir{ti} = [fish_dir 'td072817 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/NO RAD4'];
end
batch_count = batch_count+1;

B0814_RAD1_baseDir = cell(1,B0814R_max_tp);
for ti = 1:B0814R_max_tp
    B0814_RAD1_baseDir{ti} = [fish_dir 'td081417 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD1'];
end
batch_count = batch_count+1;

B0814_RAD2_baseDir = cell(1,B0814R_max_tp);
for ti = 1:B0814R_max_tp
    B0814_RAD2_baseDir{ti} = [fish_dir 'td081417 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD2'];
end
batch_count = batch_count+1;

B0814_RAD3_baseDir = cell(1,B0814R_max_tp);
for ti = 1:B0814R_max_tp
    B0814_RAD3_baseDir{ti} = [fish_dir 'td081417 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD3'];
end
batch_count = batch_count+1;

B0814_RAD4_baseDir = cell(1,B0814R_max_tp);
for ti = 1:B0814R_max_tp
    B0814_RAD4_baseDir{ti} = [fish_dir 'td081417 Adult Transplant ZMEL1 ' num2str(2*ti-1) 'DPT/RAD4'];
end
batch_count = batch_count+1;

B0903_times = [1 3 4 5 7 9 11 13];

B0903_NORAD1_baseDir = cell(1,B0903NR_max_tp);
for ti = 1:B0903NR_max_tp
    B0903_NORAD1_baseDir{ti} = [fish_dir 'td090317 Adult Transplant ZMEL1 ' num2str(B0903_times(ti)) 'DPT/NO RAD1'];
end
batch_count = batch_count+1;

B0903_NORAD2_baseDir = cell(1,B0903NR_max_tp);
for ti = 1:B0903NR_max_tp
    B0903_NORAD2_baseDir{ti} = [fish_dir 'td090317 Adult Transplant ZMEL1 ' num2str(B0903_times(ti)) 'DPT/NO RAD2'];
end
batch_count = batch_count+1;

B0903_NORAD3_baseDir = cell(1,B0903NR_max_tp);
for ti = 1:B0903NR_max_tp
    B0903_NORAD3_baseDir{ti} = [fish_dir 'td090317 Adult Transplant ZMEL1 ' num2str(B0903_times(ti)) 'DPT/NO RAD3'];
end
batch_count = batch_count+1;

B0903_NORAD4_baseDir = cell(1,B0903NR_max_tp);
for ti = 1:B0903NR_max_tp
    B0903_NORAD4_baseDir{ti} = [fish_dir 'td090317 Adult Transplant ZMEL1 ' num2str(B0903_times(ti)) 'DPT/NO RAD4'];
end
batch_count = batch_count+1;

B090619_baseDir = cell(1,B090619_max_tp);
for ti = 1:B090619_max_tp
    B090619_baseDir{ti} = [fish_dir 'td090619_' num2str(B090619_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B091119_baseDir = cell(1,B091119_max_tp);
for ti = 1:B091119_max_tp
    B091119_baseDir{ti} = [fish_dir 'td091119_' num2str(B091119_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B091919_baseDir = cell(1,B091919_max_tp);
for ti = 1:B091919_max_tp
    B091919_baseDir{ti} = [fish_dir 'td091919_' num2str(B091919_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B092319_baseDir = cell(1,B092319_max_tp);
for ti = 1:B092319_max_tp
    B092319_baseDir{ti} = [fish_dir 'td092319_' num2str(B092319_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B092719_baseDir = cell(1,B092719_max_tp);
for ti = 1:B092719_max_tp
    B092719_baseDir{ti} = [fish_dir 'td092719_' num2str(B092719_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B100119_baseDir = cell(1,B100119_max_tp);
for ti = 1:B100119_max_tp
    B100119_baseDir{ti} = [fish_dir 'td100119_' num2str(B100119_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B100919_baseDir = cell(1,B100919_max_tp);
for ti = 1:B100919_max_tp
    B100919_baseDir{ti} = [fish_dir 'td100919_' num2str(B100919_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B101319_baseDir = cell(1,B101319_max_tp);
for ti = 1:B101319_max_tp
    B101319_baseDir{ti} = [fish_dir 'td101319_' num2str(B101319_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B101719_baseDir = cell(1,B101719_max_tp);
for ti = 1:B101719_max_tp
    B101719_baseDir{ti} = [fish_dir 'td101719_' num2str(B101719_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B102119_baseDir = cell(1,B102119_max_tp);
for ti = 1:B102119_max_tp
    B102119_baseDir{ti} = [fish_dir 'td102119_' num2str(B102119_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

B102519_baseDir = cell(1,B102519_max_tp);
for ti = 1:B102519_max_tp
    B102519_baseDir{ti} = [fish_dir 'td102519_' num2str(B102519_times(ti)) 'DPT/NO RAD'];
end
batch_count = batch_count+1;

% gender files
B0130_RAD_gender = [fish_dir 'gender_td013017_RAD.txt'];
B0130_NORAD_gender = [fish_dir 'gender_td013017_NORAD.txt'];
B0213_RAD_gender = [fish_dir 'gender_td021317_RAD.txt'];
B0213_NORAD_gender = [fish_dir 'gender_td021317_NORAD.txt'];
B0709_RAD1_gender = [fish_dir 'gender_td070917_RAD1.txt'];
B0709_RAD2_gender = [fish_dir 'gender_td070917_RAD2.txt'];
B0728_NORAD1_gender = [fish_dir 'gender_td072817_NORAD1.txt'];
B0728_NORAD2_gender = [fish_dir 'gender_td072817_NORAD2.txt'];
B0728_NORAD3_gender = [fish_dir 'gender_td072817_NORAD3.txt'];
B0728_NORAD4_gender = [fish_dir 'gender_td072817_NORAD4.txt'];
B0814_RAD1_gender = [fish_dir 'gender_td081417_RAD1.txt'];
B0814_RAD2_gender = [fish_dir 'gender_td081417_RAD2.txt'];
B0814_RAD3_gender = [fish_dir 'gender_td081417_RAD3.txt'];
B0814_RAD4_gender = [fish_dir 'gender_td081417_RAD4.txt'];
B0903_NORAD1_gender = [fish_dir 'gender_td090317_NORAD1.txt'];
B0903_NORAD2_gender = [fish_dir 'gender_td090317_NORAD2.txt'];
B0903_NORAD3_gender = [fish_dir 'gender_td090317_NORAD3.txt'];
B0903_NORAD4_gender = [fish_dir 'gender_td090317_NORAD4.txt'];
B090619_NORAD_gender = [fish_dir 'gender_td090619_NORAD.txt'];
B091119_NORAD_gender = [fish_dir 'gender_td091119_NORAD.txt'];
B091919_NORAD_gender = [fish_dir 'gender_td091919_NORAD.txt'];
B092319_NORAD_gender = [fish_dir 'gender_td092319_NORAD.txt'];
B092719_NORAD_gender = [fish_dir 'gender_td092719_NORAD.txt'];
B100119_NORAD_gender = [fish_dir 'gender_td100119_NORAD.txt'];
B100919_NORAD_gender = [fish_dir 'gender_td100919_NORAD.txt'];
B101319_NORAD_gender = [fish_dir 'gender_td101319_NORAD.txt'];
B101719_NORAD_gender = [fish_dir 'gender_td101719_NORAD.txt'];
B102119_NORAD_gender = [fish_dir 'gender_td102119_NORAD.txt'];
B102519_NORAD_gender = [fish_dir 'gender_td102519_NORAD.txt'];

% lookup files
B0130_RAD_lookup = [fish_dir 'lookup_td013017_RAD.txt'];
B0130_NORAD_lookup = [fish_dir 'lookup_td013017_NORAD.txt'];
B0213_RAD_lookup = [fish_dir 'lookup_td021317_RAD.txt'];
B0213_NORAD_lookup = [fish_dir 'lookup_td021317_NORAD.txt'];
B0709_RAD1_lookup = [fish_dir 'lookup_td070917_RAD1.txt'];
B0709_RAD2_lookup = [fish_dir 'lookup_td070917_RAD2.txt'];
B0728_NORAD1_lookup = [fish_dir 'lookup_td072817_NORAD1.txt'];
B0728_NORAD2_lookup = [fish_dir 'lookup_td072817_NORAD2.txt'];
B0728_NORAD3_lookup = [fish_dir 'lookup_td072817_NORAD3.txt'];
B0728_NORAD4_lookup = [fish_dir 'lookup_td072817_NORAD4.txt'];
B0814_RAD1_lookup = [fish_dir 'lookup_td081417_RAD1.txt'];
B0814_RAD2_lookup = [fish_dir 'lookup_td081417_RAD2.txt'];
B0814_RAD3_lookup = [fish_dir 'lookup_td081417_RAD3.txt'];
B0814_RAD4_lookup = [fish_dir 'lookup_td081417_RAD4.txt'];
B0903_NORAD1_lookup = [fish_dir 'lookup_td090317_NORAD1.txt'];
B0903_NORAD2_lookup = [fish_dir 'lookup_td090317_NORAD2.txt'];
B0903_NORAD3_lookup = [fish_dir 'lookup_td090317_NORAD3.txt'];
B0903_NORAD4_lookup = [fish_dir 'lookup_td090317_NORAD4.txt'];
B090619_NORAD_lookup = [fish_dir 'lookup_td090619_NORAD.txt'];
B091119_NORAD_lookup = [fish_dir 'lookup_td091119_NORAD.txt'];
B091919_NORAD_lookup = [fish_dir 'lookup_td091919_NORAD.txt'];
% B092319_NORAD_lookup;
B092719_NORAD_lookup = [fish_dir 'lookup_td092719_NORAD.txt'];
B100119_NORAD_lookup = [fish_dir 'lookup_td100119_NORAD.txt'];
B100919_NORAD_lookup = [fish_dir 'lookup_td100919_NORAD.txt'];
B101319_NORAD_lookup = [fish_dir 'lookup_td101319_NORAD.txt'];
B101719_NORAD_lookup = [fish_dir 'lookup_td101719_NORAD.txt'];
B102119_NORAD_lookup = [fish_dir 'lookup_td102119_NORAD.txt'];
B102519_NORAD_lookup = [fish_dir 'lookup_td102519_NORAD.txt'];

rad_batch = [1 3 5:6 11:14];
norad_batch = [2 4 7:10 15:29];

% will iterate through these states
batches = {'B013017R1','B013017NR1','B021317R1','B021317NR1','B070917R1','B070917R2',...
    'B072817NR1','B072817NR2','B072817NR3','B072817NR4','B081417R1',...
    'B081417R2','B081417R3','B081417R4','B090317NR1',...
    'B090317NR2','B090317NR3','B090317NR4','B090619NR','B091119NR','B091919NR',...
    'B092319NR','B092719NR','B100119NR','B100919NR','B101319NR','B101719NR',...
    'B102119NR','B102519NR'}; 

batch_dir = {B0130_RAD_baseDir,B0130_NORAD_baseDir,B0213_RAD_baseDir,B0213_NORAD_baseDir,...
    B0709_RAD_baseDir,B0709_RAD2_baseDir,B0728_NORAD1_baseDir,...
    B0728_NORAD2_baseDir,B0728_NORAD3_baseDir,B0728_NORAD4_baseDir,...
    B0814_RAD1_baseDir,B0814_RAD2_baseDir,B0814_RAD3_baseDir,...
    B0814_RAD4_baseDir,B0903_NORAD1_baseDir,B0903_NORAD2_baseDir,...
    B0903_NORAD3_baseDir,B0903_NORAD4_baseDir,B090619_baseDir,...
    B091119_baseDir,B091919_baseDir,B092319_baseDir,B092719_baseDir,...
    B100119_baseDir,B100919_baseDir,B101319_baseDir,B101719_baseDir,...
    B102119_baseDir,B102519_baseDir};

exposure_times = {B0130_exp_times,B0130_exp_times,B0213_exp_times,B0213_exp_times,...
    B0709_exp_times,B0709_exp_times,B0728_exp_times,B0728_exp_times,...
    B0728_exp_times,B0728_exp_times,B0814_exp_times,B0814_exp_times,...
    B0814_exp_times,B0814_exp_times,B0903_exp_times,B0903_exp_times,...
    B0903_exp_times,B0903_exp_times,B090619_exp_times,B091119_exp_times,...
    B091919_exp_times,B092319_exp_times,B092719_exp_times,B100119_exp_times,...
    B100919_exp_times,B101319_exp_times,B101719_exp_times,B102119_exp_times,...
    B102519_exp_times};

% naming_conventions = {B0130_img_naming_conv,B0130_img_naming_conv,B0130_img_naming_conv,...
%     B0130_img_naming_conv,B0709_img_naming_conv,B0709_img_naming_conv,...
%     B0709_img_naming_conv,B0709_img_naming_conv,B0709_img_naming_conv,...
%     B0709_img_naming_conv,B0709_img_naming_conv,B0709_img_naming_conv,...
%     B0709_img_naming_conv,B0709_img_naming_conv,B0709_img_naming_conv,...
%     B0709_img_naming_conv,B0709_img_naming_conv,B0709_img_naming_conv,...
%     B0709_img_naming_conv,B0709_img_naming_conv,B0709_img_naming_conv,...
%     B0709_img_naming_conv};
naming_conventions = cell(1,batch_count);
for bc = 1:4
    naming_conventions{bc} = B0130_img_naming_conv;
end
for bc = 5:batch_count
    naming_conventions{bc} = B0709_img_naming_conv;
end

gender_files = {B0130_RAD_gender,B0130_NORAD_gender,B0213_RAD_gender,B0213_NORAD_gender,...
    B0709_RAD1_gender,B0709_RAD2_gender,B0728_NORAD1_gender,...
    B0728_NORAD2_gender,B0728_NORAD3_gender,B0728_NORAD4_gender,...
    B0814_RAD1_gender,B0814_RAD2_gender,B0814_RAD3_gender,...
    B0814_RAD4_gender,B0903_NORAD1_gender,B0903_NORAD2_gender,...
    B0903_NORAD3_gender,B0903_NORAD4_gender,B090619_NORAD_gender,...
    B091119_NORAD_gender,B091919_NORAD_gender,B092319_NORAD_gender,...
    B092719_NORAD_gender,B100119_NORAD_gender,B100919_NORAD_gender,...
    B101319_NORAD_gender,B101719_NORAD_gender,B102119_NORAD_gender,...
    B102519_NORAD_gender};

lookup_files = {B0130_RAD_lookup,B0130_NORAD_lookup,B0213_RAD_lookup,B0213_NORAD_lookup,...
    B0709_RAD1_lookup,B0709_RAD2_lookup,B0728_NORAD1_lookup,...
    B0728_NORAD2_lookup,B0728_NORAD3_lookup,B0728_NORAD4_lookup,...
    B0814_RAD1_lookup,B0814_RAD2_lookup,B0814_RAD3_lookup,...
    B0814_RAD4_lookup,B0903_NORAD1_lookup,B0903_NORAD2_lookup,...
    B0903_NORAD3_lookup,B0903_NORAD4_lookup,B090619_NORAD_lookup,'',...
    B091919_NORAD_lookup,'',B092719_NORAD_lookup,B100119_NORAD_lookup,...
    B100919_NORAD_lookup,B101319_NORAD_lookup,B101719_NORAD_lookup,...
    B102119_NORAD_lookup,B102519_NORAD_lookup};

elimination_files = cell(1,batch_count);
for bi = 1:batch_count
    elimination_files{bi} = '';
end

implantSizes = {'5x10^5','5x10^5','5x10^5','5x10^5','1x10^5','1x10^5','1x10^5',...
    '1x10^5','1x10^5','1x10^5','1x10^5','1x10^5','1x10^5','1x10^5','1x10^5',...
    '1x10^5','1x10^5','1x10^5','1x10^6','1x10^6','5x10^6','5x10^6','5x10^6',...
    '5x10^6','5x10^6','1x10^7','5x10^6','1x10^7','1x10^7'}; % first half of group name folder names
implantSizeNames = {'5e5','5e5','5e5','5e5','1e5','1e5','1e5','1e5','1e5',...
    '1e5','1e5','1e5','1e5','1e5','1e5','1e5','1e5','1e5','1e6','1e6','5e6',...
    '5e6','5e6','5e6','5e6','1e7','5e6','1e7','1e7'};

locations = cell(1,batch_count);
for bi = 1:25
    locations{bi} = 'VENTRAL';
end
locations{26} = 'MIXED';
locations{27} = 'DORSAL';
locations{28} = 'MIXED';
locations{29} = 'MIXED';

timePoints = cell(1,batch_count);
for bi = 1:14
    timePoints{bi} = complete_times(1:length(batch_dir{bi}));
end
% last for have a time point at day 4
for bi = 15:18
    timePoints{bi} = [complete_times(1:2) '4' complete_times(3:length(batch_dir{bi})-1)];
end
timePoints{19} = num2cell(B090619_times);
timePoints{20} = num2cell(B091119_times);
timePoints{21} = num2cell(B091919_times);
timePoints{22} = num2cell(B092319_times);
timePoints{23} = num2cell(B092719_times);
timePoints{24} = num2cell(B100119_times);
timePoints{25} = num2cell(B100919_times);
timePoints{26} = num2cell(B101319_times);
timePoints{27} = num2cell(B101719_times);
timePoints{28} = num2cell(B102119_times);
timePoints{29} = num2cell(B102519_times);

def_nose = [40 300];
def_eye = def_nose + [90 -10];
fish_specs_default = struct('minimum_area',50000,'maximum_length',1380,...
    'maximum_width',700,'nose_coordinate',def_nose,'eye_diameter',81,...
    'eye_coordinate',def_eye,'clustering_threshold',10,...
    'experimental_exposure_times',B0130_exp_times,'imageConversionFactor',257,...
    'eyeBrightnessFactor',30,'rfp_gfp_ratio',2,'minimum_eye_area',2500,...
    'mean_intensities',[6.22e4,2.00e3,2.00e3,5.13e4],'bw_threshold_high',0.04,...
    'bw_threshold_low',0.015,'normal_exposure_times',normal_exposure_times,...
    'reference_threshold',0.4,'clustering_threshold_day1',10,...
    'edge_threshold',0.0012,'min_pigmentation_threshold', 0.01,...
    'pigmentation_use_index',2,'selector',[]);

% setting to display results of image analysis
start_step = 1;
fish_to_skip = [];
start_batch = 1;
show = 0;

% find finished runs
finished_tumor_files = dir([baseName '*summary_tumors.mat']);
finishedIDs = cell(1,length(finished_tumor_files)); 
for ftf = 1:length(finished_tumor_files)
    finishedIDs{ftf} = finished_tumor_files(ftf).name(1:end-19);
end

selectors = cell(1,length(batches));

part1 = norad_batch(1:5);
part2 = norad_batch(6:10);
part3 = rad_batch(1:4);
part4 = rad_batch(5:8);
part5 = norad_batch([11:13 15]);
part6 = norad_batch(16:19);
part7 = norad_batch([14 20:21]);

run_list = [rad_batch norad_batch];
pf_specs = cell(1,length(run_list));
for ss = 1:length(pf_specs)
    new_specs = fish_specs_default;
    new_specs.experimental_exposure_times = exposure_times{run_list(ss)};
    new_specs.selector = selectors{run_list(ss)};
    pf_specs{ss} = new_specs;
end

pf_batches = batches(run_list);
pf_batch_dir = batch_dir(run_list);
pf_timePoints = timePoints(run_list);
pf_naming_conventions = naming_conventions(run_list);
pf_gender_files = gender_files(run_list);
pf_lookup_files = lookup_files(run_list);
pf_elimination_files = elimination_files(run_list);
pf_implantSizes = implantSizes(run_list);
pf_implantSizeNames = implantSizeNames(run_list);
pf_locations = locations(run_list);

% save results
tumor_lists = cell(1,length(run_list));
total_lists = cell(1,length(run_list));
parfor ss = 1:length(run_list)
    % set up variables
    batch = pf_batches(ss);
    baseDir = pf_batch_dir{ss};
    timePoint = pf_timePoints{ss};
    fish_specs = pf_specs{ss};
    img_names = pf_naming_conventions{ss};
    gender_file = pf_gender_files{ss};
    lookup_file = pf_lookup_files{ss};
    elimination_file = pf_elimination_files{ss};
    implantSize = pf_implantSizes(ss);
    implantSizeName = pf_implantSizeNames(ss);
    location = pf_locations(ss);
    fish_struct_save_name = [baseName '_' batch{1} '_' ...
        location{1} '_' implantSizeName{1}];
    disp('--->')
    if(any(ismember(finishedIDs,fish_struct_save_name)))
        disp(['Skipping ' fish_struct_save_name '...']);
    else
        % clc;
        disp(['Running ' fish_struct_save_name '...']);
        % run code
        [tumor_lists{ss},total_lists{ss}] = CreateFishStructNested_v4_1(start_step,fish_to_skip,show,batch,baseDir,timePoint,fish_specs,...
            img_names,gender_file,lookup_file,elimination_file,implantSize,implantSizeName,location,fish_struct_save_name);
    end
end
            