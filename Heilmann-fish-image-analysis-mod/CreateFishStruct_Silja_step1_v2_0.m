%% CreateFishStruct_Silja_step1_v2_0
%  Version 2.0
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 3/12/18
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Version History
%  2.0: using red rgb image as bf and storing rgb as rgb separately

function [BYTIME] = CreateFishStruct_Silja_step1_v2_0(baseDir,timePoint,batch,implantSize,implantSizeName,location)
%% Add original images to structure array "FISH"
%  As the data structure FISH is created each fish get a unique ID contaning
%  information about  original fish number, batch, implant size and implant
%  location. Original images are then imported from their folders into matlab and
%  saved in FISH 

tic;

disp('Step 1: Initializing the data structure FISH and reading in original images...');

%Initiale data structure
BYTIME(length(timePoint)).fishes = struct();
nfsh_init = 200;
FISH(nfsh_init).fishID = 'fID';
counter = zeros(length(timePoint),nfsh_init);
max_fish = 0;

for tt = 1:length(timePoint)
    clear FISH;
    FISH(nfsh_init).fishID = 'fID';
    for bb = 1:length(batch)
        for ii = 1:length(implantSize)
            for ll = 1:length(location)
                % make list of fish images
                fileListFISH  = dir([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/*.tif']);
                
                % sort images by name
                % bf1, bf2, gfp1, gfp2, rfp1, rfp2
                
                for f = 1:length(fileListFISH) % iterate through fish images
                    
                    % need to account for 6 images for each fish
                    img_name = fileListFISH(f).name(1:end-4);
                    if(img_name(end-2:end-1)=='bf')
                        % name is #bf#
                        sort_name = img_name(end-2:end);
                        fnum = str2double(img_name(1:end-3));
                    elseif(img_name(end-2:end-1)=='fp')
                        % name is #gfp# or #rfp#
                        sort_name = img_name(end-3:end);
                        fnum = str2double(img_name(1:end-4));
                    else
                        % un-used image
                        sort_name = 'un-used';
                    end
                    
                    switch sort_name
                        case 'bf1'
                            fishID = ['fish_' num2str(fnum) '_' batch{bb} '_' implantSizeName{ii} '_' location{ll}];
                            fish = fnum; % these pictures should be labeled with a fish number to identify the fish
                            max_fish = max(fish,max_fish);
%                             if(isfield(lookup,fishID))
%                                 fish = lookup.(fishID);
%                             else
%                                 fish = fnum;
%                                 lookup.(fishID) = fish;
%                             end
                            % save info about new fish
                            FISH(fish).fishID       = fishID;
                            FISH(fish).batch        = batch{bb};
                            FISH(fish).implantSize  = implantSize{ii};
                            FISH(fish).location     = location{ll};
                            FISH(fish).time = timePoint{tt};
                            counter(tt,fish) = counter(tt,fish)+1;
                            FISH(fish).rgb1 = imread([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name]);
                            FISH(fish).bf1 = FISH(fish).rgb1;
                        case 'bf2'
                            FISH(fish).rgb2 = imread([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name]);
                            FISH(fish).bf2 = FISH(fish).rgb2;
                        case 'gfp1'
                            FISH(fish).gfp1 = imread([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name]);
                        case 'gfp2'
                            FISH(fish).gfp2 = imread([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name]);
                        case 'rfp1'
                            FISH(fish).rfp1 = imread([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name]);
                        case 'rfp2'
                            FISH(fish).rfp2 = imread([baseDir '/'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name]);
                        case 'un-used'
                            disp(['Image named "' img_name '" was discarded from /'  batch{bb}  '-renamed/' timePoint{tt} '/' implantSize{ii} '' location{ll}]);
                    end %
                end
                
            end
        end
    end
    BYTIME(tt).fishes = FISH(1:max_fish);
end

% find fish with pictures at all times
if(length(timePoint)>1)
    dataAllTimes = find(sum(counter)==length(timePoint));
else
    dataAllTimes = find(counter);
end
disp([num2str(length(dataAllTimes)) ' fish out of ' num2str(length(BYTIME(1).fishes)) ' with all data found']);

% list fish with data at specific times
dataAt = cell(1,length(timePoint));
for i=1:length(timePoint)
   dataAt{i} = find(counter(i,:));
end

BYTIME(end).dataAllTimes = dataAllTimes;
BYTIME(end).dataAt = dataAt;

% this links the labels of the fish at the last time point to labels at
% previous time points (not useful here since fish identities are known).
% This is here for compatibility to prevent errors later
onevec = ones(1,length(BYTIME));
for ff=1:length(BYTIME(end).fishes)
    BYTIME(end).fishes(ff).config = ff*onevec;
end

toc
disp(' ')