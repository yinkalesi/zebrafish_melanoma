%% reviseTumorImageData_v1_2
%  version 1.2
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 12/09/19

%% Version History
%  1.1: using getL_Revised_v1_2
%  1.2: using getL_Revised_v1_3 -> using revised files to facilitate making
%  minor corrections
%  (2/20/20): switched to getL_Revised_v1_4 (changes picture displayed)

function [rev_tum,rev_Ldivs] = reviseTumorImageData_v1_2(tum_file,img_file,new_name,old_rev)

% load files
temp = load(tum_file);
tumors = temp.tumors;
temp = load(img_file);
images = temp.images;

% get user input
if(exist('old_rev','var'))
    old_rev_dat = load(old_rev);
    [rev_tum,rev_Ldivs] = getL_Revised_v1_4(tumors,images,old_rev_dat.rev_tum);
else
    [rev_tum,rev_Ldivs] = getL_Revised_v1_4(tumors,images);
end

% save files
save([new_name '_tumors_revised.mat'],'rev_tum');
save([new_name '_Ldivs_revised.mat'],'rev_Ldivs');

end

