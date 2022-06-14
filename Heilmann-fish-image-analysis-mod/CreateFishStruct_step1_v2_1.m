%% CreateFishStruct_step1
%  Version 2.1
%  Author: Dr. Silja Heilmann
%  Contributors: Adeyinka Lesi
%  Date: 2/09/20
%  Project: Tumor Growth, Logarithmic Continuum Form
% Add original images to structure array "FISH"
%  As the data structure FISH is created each fish get a unique ID containing
%  information about  original fish number, batch, implant size and implant
%  location. Original images are then imported from their folders into MATLAB and
%  saved in FISH 
%% Version History
%  2.1: adding selector to pick which fish images are analyzed

function [BYTIME] = CreateFishStruct_step1_v2_1(baseDir,timePoint,batch,implantSize,location,img_names,specs)

tic;

disp('Step 1: Initializing the data structure FISH and reading in original images...');

fileListFISH = dir([baseDir{1} '/*tif']);
% assuming 10 images per fish
nfsh_init = round(0.1*length(fileListFISH));

if(~isempty(specs.selector))
	selector = specs.selector;
    time_list = 1:length(selector);
else
	selector = [];
    time_list = 1:length(timePoint);
end

% Initialize data structure
BYTIME(length(time_list)).fishes = struct();
FISH(nfsh_init).fishID = zeros(1,nfsh_init);
uniid = 0; % counter for counting fish as they are added to data structure


for bb = 1:length(batch)
    for ii = 1:length(implantSize)
        for ll = 1:length(location)
            for tt = time_list
        
                % make list of fish images
                fileListFISH  = dir([baseDir{tt} '/*tif']);
                if(isempty(fileListFISH))
                    % B013017 15DPT is made of jpg images
                    fileListFISH  = dir([baseDir{tt} '/*jpg']);
                end
                if(isempty(fileListFISH))
                    error(['Cannot find any files in directory: ',...
                        baseDir{tt}]);
                end
                
                % need to figure out true order of files - the files are
                % ordered by number, but matlab orders it so that '100' is
                % ahead of '65' - so I will extract the numbers from the
                % names and sort
                fileNums = zeros(1,length(fileListFISH));
                for f = 1:length(fileListFISH)
                    fname = fileListFISH(f).name;
                    % find underscore position
                    usi = find(fname=='_',1);
                    % find hyphen position
                    hi = find(fname=='-',1);
                    if(isempty(usi))
                        usi = find(fname=='.',1); % one of the color files with no underscore
                    end
                    fileNums(f) = str2double(fname(hi+1:usi-1));
                end
                [~,fileOrder] = sort(fileNums);
                
                % rubric for naming
                % Color: Snap-#.tif
                clrid = img_names.color;
                % bf: Snap-#_TL Brightfield_1.tif
                bfid = img_names.brightfield;
                % gfp: Snap-#_EGFP_3.tif
                gfpid = img_names.gfp;
                % rfp: Snap-#_Rhodamine_2.tif
                rfpid = img_names.rfp;
                % overlay: Snap-#_c1-3.tif
                ovlid = img_names.overlay;
                
                % sort by name;
                img_num = 0;
                thresh = 10;
                for flf = 1:length(fileListFISH) % iterate through fish images
                        f = fileOrder(flf);
                        fname = fileListFISH(f).name;
                        % find underscore position
                        usi = find(fname=='_',1);
                        % find hyphen position
                        hi = find(fname=='-',1);
                        % period position
                        pi = find(fname=='.',1);
                        if(isempty(usi))
                            usi = pi; % one of the color files with no underscore
                        end
                        imgid = fname(usi:pi);
                        img_num = mod(img_num,thresh);
%                         disp(['(' num2str(f) ',' num2str(uniid) '): ' imgid '; ' num2str(img_num)]);
                        if(img_num == 0)
                            % new fish
                            uniid = uniid + 1;
                            fishID = ['fish' num2str(uniid) '_t' num2str(tt)...
                                '_' batch{bb} '_' implantSize{ii} '_' location{ll}];
                            % save info about new fish
                            FISH(uniid).fishID       = fishID;
                            FISH(uniid).batch        = batch{bb};
                            FISH(uniid).implantSize  = implantSize{ii};
                            FISH(uniid).location     = location{ll};
                            FISH(uniid).time = timePoint{tt};
                            FISH(uniid).snap = zeros(1,thresh);
                        end
                        switch imgid
                            case clrid
                                img_num = img_num + 1;
                                if(img_num > thresh*0.5)
                                    FISH(uniid).rgb1 = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                else
                                    FISH(uniid).rgb2 = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                end
                                FISH(uniid).snap(img_num) = str2double(fname(hi+1:usi-1));
                            case bfid
                                img_num = img_num + 1;
                                if(img_num > thresh*0.5)
                                    temp = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                    FISH(uniid).bf1 = temp(:,:,1);
                                else
                                    temp = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                    FISH(uniid).bf2 = temp(:,:,1);
                                end
                                FISH(uniid).snap(img_num) = str2double(fname(hi+1:usi-1));
                            case gfpid
                                img_num = img_num + 1;
                                if(img_num > thresh*0.5)
                                    temp = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                    FISH(uniid).gfp1 = temp(:,:,2);
                                else
                                    temp = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                    FISH(uniid).gfp2 = temp(:,:,2);
                                end
                                FISH(uniid).snap(img_num) = str2double(fname(hi+1:usi-1));
                            case rfpid
                                img_num = img_num + 1;
                                if(img_num > thresh*0.5)
                                    temp = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                    FISH(uniid).rfp1 = temp(:,:,1);
                                else
                                    temp = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                    FISH(uniid).rfp2 = temp(:,:,1);
                                end
                                FISH(uniid).snap(img_num) = str2double(fname(hi+1:usi-1));
                            case ovlid
                                % should be an rgb image; check
                                img_num = img_num + 1;
                                if(img_num > thresh*0.5)
                                    FISH(uniid).ovl1 = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                else
                                    FISH(uniid).ovl2 = imread([baseDir{tt} '/' fileListFISH(f).name]);
                                end
                                FISH(uniid).snap(img_num) = str2double(fname(hi+1:usi-1));
                            otherwise
                                % this file is not viable data
                                if(img_num == 0)
                                    % fix id situation if this is meant to
                                    % be first image
                                    uniid = uniid - 1;
                                end
                        end %
                end
				if(~isempty(selector))
					BYTIME(tt).fishes = FISH(selector{tt});
				else
					BYTIME(tt).fishes = FISH(1:uniid);
				end
                uniid = 0;
            end
        end
    end
end

dataAt = cell(1,length(time_list));
for tt = 1:length(dataAt)
    dataAt{tt} = 1:length(BYTIME(tt).fishes);
end
BYTIME(end).dataAt = dataAt;

toc
disp(' ')