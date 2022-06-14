%% runCreateFishStruct_Silja_v2_0
%  Version 2.0
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 2/21/20
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Version History
%  1.1: changing day1 clustering/alignment issues/better coverage of
%  pigmented tumors
%  2.0: updating (2/21/20)

% clear variables;
baseName = 'Silja_022520_r2';
fish_dir = 'C:\Users\David\Documents\Yinka\fish_images\';

timePoint_standard = {'Day1','Week1','Week2'};

% exposure times ([BF;GFP;RFP;COLOR])
exp_times = [2 2 2;...
    300 300 300;...
    300 300 300;...
    2 2 2];

normal_exposure_times = [2; 300; 300; 2];

% picture naming convention
img_names_standard = struct('color','.','brightfield','.',...
    'gfp','.','rfp','.','overlay','.')';

% will iterate through these states
batches = {'Batch1','Batch1','Batch1','Batch1','Batch1','Batch1',...
    'Batch2','Batch2','Batch2','Batch2','Batch2','Batch2'};
% batches   = {'B1_IM1_V','B1_IM2_V','B1_IM3_V','B1_IM1_D','B1_IM2_D',...
%     'B1_IM3_D','B2_IM1_V','B2_IM2_V','B2_IM3_V','B2_IM1_D','B2_IM2_D','B2_IM3_D'}; %names of batch folders
batch_count = 12;
implantSizes = {'1x10^5','5x10^5','1x10^6','1x10^5','5x10^5','1x10^6',...
    '1x10^5','5x10^5','1x10^6','1x10^5','5x10^5','1x10^6'}; % first half of group name folder names
implantSizeNames = {'1e5','5e5','1e6','1e5','5e5','1e6',...
    '1e5','5e5','1e6','1e5','5e5','1e6'};
locations = {'VENTRAL','VENTRAL','VENTRAL','DORSAL','DORSAL','DORSAL',...
    'VENTRAL','VENTRAL','VENTRAL','DORSAL','DORSAL','DORSAL'}; % second half of group name folder names

batch_dir = cell(1,batch_count);
for bi = 1:batch_count
    batch_dir{bi} = fish_dir;
end

timePoints = cell(1,batch_count);
for bi = 1:batch_count
    timePoints{bi} = timePoint_standard;
end

exposure_times = cell(1,batch_count);
for bi = 1:batch_count
    exposure_times{bi} = exp_times;
end

naming_conventions = cell(1,batch_count);
for bi = 1:batch_count
    naming_conventions{bi} = img_names_standard;
end

% gender files
batch1_im1_ventral = [fish_dir 'gender.txt'];
batch1_im2_ventral = [fish_dir 'gender.txt'];
batch1_im3_ventral = [fish_dir 'gender.txt'];
batch1_im1_dorsal = [fish_dir 'gender.txt'];
batch1_im2_dorsal = [fish_dir 'gender.txt'];
batch1_im3_dorsal = [fish_dir 'gender.txt'];
batch2_im1_ventral = [fish_dir 'gender.txt'];
batch2_im2_ventral = [fish_dir 'gender.txt'];
batch2_im3_ventral = [fish_dir 'gender.txt'];
batch2_im1_dorsal = [fish_dir 'gender.txt'];
batch2_im2_dorsal = [fish_dir 'gender.txt'];
batch2_im3_dorsal = [fish_dir 'gender.txt'];

gender_files = {batch1_im1_ventral,batch1_im2_ventral,batch1_im3_ventral,...
    batch1_im1_dorsal,batch1_im2_dorsal,batch1_im3_dorsal,...
    batch2_im1_ventral,batch2_im2_ventral,batch2_im3_ventral,...
    batch2_im1_dorsal,batch2_im2_dorsal,batch2_im3_dorsal};

elimination_files = cell(1,batch_count);
for bi = 1:batch_count
    elimination_files{bi} = '';
end

def_nose = [10 175];
def_eye = def_nose + [60 -5];
fish_specs_default = struct('minimum_area',35000,'maximum_length',1000,...
    'maximum_width',350,'nose_coordinate',def_nose,'eye_diameter',54,...
    'eye_coordinate',def_eye,'clustering_threshold',10,...
    'experimental_exposure_times',exp_times,'imageConversionFactor',1,...
    'eyeBrightnessFactor',10,'rfp_gfp_ratio',2,'minimum_eye_area',1500,...
    'mean_intensities',[6.22e4,2.00e3,2.00e3],'bw_threshold_high',0.05,...
    'bw_threshold_low',0.01,'normal_exposure_times',normal_exposure_times,...
    'reference_threshold',0.2,'clustering_threshold_day1',10,...
    'edge_threshold',0.015,'min_pigmentation_threshold', 0.01,...
    'pigmentation_use_index',2,'selector',[]);

% setting to display results of image analysis
start_step = 1;
fish_to_skip = [];
start_batch = 1;
show = 0;

% finished runs
known_summary_files = dir([baseName '_*summary_tumors.mat']);
finishedIDs = cell(1,length(known_summary_files));
for i = 1:length(known_summary_files)
    finishedIDs{i} = known_summary_files(i).name(1:find(known_summary_files(i).name=='.',1,'last')-16);
end

selectors = cell(1,length(batches));

run_list = 1:batch_count;
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
% pf_lookup_files = lookup_files(run_list);
pf_elimination_files = elimination_files(run_list);
pf_implantSizes = implantSizes(run_list);
pf_implantSizeNames = implantSizeNames(run_list);
pf_locations = locations(run_list);



% save results
tumor_lists = cell(1,length(run_list));
total_lists = cell(1,length(run_list));
for ss = 1:length(run_list)
    % set up variables
    batch = pf_batches(ss);
    baseDir = pf_batch_dir{ss};
    timePoint = pf_timePoints{ss};
    fish_specs = pf_specs{ss};
    img_names = pf_naming_conventions{ss};
    gender_file = pf_gender_files{ss};
    %     lookup_file = pf_lookup_files{ss};
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
%         clc;
        disp(['Running ' fish_struct_save_name '...']);
        % run code
        [tumor_lists{ss},total_lists{ss}] = CreateFishStructNested_Silja_v2_0(start_step,fish_to_skip,show,batch,baseDir,timePoint,fish_specs,...
            img_names,gender_file,elimination_file,implantSize,implantSizeName,location,fish_struct_save_name);
    end
end

            