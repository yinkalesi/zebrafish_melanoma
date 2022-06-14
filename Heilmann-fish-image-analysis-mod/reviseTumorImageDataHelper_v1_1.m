%% reviseTumorImageDataHelper_v1_1
%  version 1.1
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 12/9/19
%% Version History
%  1.0: 2017 data
%  1.1: made revisions easier using reviseTumorImageData_v1_2

dat_dir = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\Continuous Tumor Growth Model\Heilmann-fish-image-analysis-mod\'; 
summary_files = dir([dat_dir 'update1101_newseg_*17*e5*tumors.mat']);
image_files =  dir([dat_dir 'update1101_newseg_*17*e5*images.mat']);

for i = 1:length(summary_files)
    tum_file = [dat_dir summary_files(i).name];
    img_file = [dat_dir image_files(i).name];
    endloc = find(tum_file=='_',1,'last')-1;
    newname = tum_file(1:endloc);
    slashloc = find(tum_file=='\',1,'last');
    if(slashloc)
        shortname = tum_file(slashloc+1:end);
    else
        shortname = tum_file;
    end
    fprintf('%i: %s\n',i,shortname);
    reviseTumorImageData_v1_2(tum_file,img_file,newname);
    move_on = questdlg(['Move on with ' shortname '?']);
    switch move_on
        case 'Yes'
            fprintf('Finished %s\n',i,shortname);
        case 'No'
            fprintf('Redoing %s\n',i,shortname);
            reviseTumorImageData_v1_2(tum_file,img_file,newname,[new_name '_tumors_revised.mat']);
            fprintf('Finished %s\n',i,shortname);
    end
end

