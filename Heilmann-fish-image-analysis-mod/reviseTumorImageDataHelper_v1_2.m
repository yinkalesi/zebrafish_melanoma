%% reviseTumorImageDataHelper_v1_2
%  version 1.2
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 1/30/20
%% Version History
%  1.0: 2017 data
%  1.1: made revisions easier using reviseTumorImageData_v1_2
%  1.2: cluster threshold back to 10

dat_dir = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\Continuous Tumor Growth Model\Heilmann-fish-image-analysis-mod\'; 
summary_files = dir([dat_dir 'ct10_lt015_ht04_021720_nosat_*19NR*tumors.mat']);
image_files =  dir([dat_dir 'ct10_lt015_ht04_021720_nosat_*19NR*images.mat']);

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
    if(~exist([newname '_tumors_revised_test1.mat'],'file'))
        reviseTumorImageData_v1_2(tum_file,img_file,newname);
    else
        fprintf('File %s already exists\n',[newname '_tumors_revised_test1.mat']);
    end
    move_on = questdlg(['Move on with ' shortname '?']);
    switch move_on
        case 'Yes'
            fprintf('Finished %i: %s\n',i,shortname);
        case 'No'
            fprintf('Redoing %i: %s\n',i,shortname);
            reviseTumorImageData_v1_2(tum_file,img_file,newname,[newname '_tumors_revised_test1.mat']);
            fprintf('Finished %i: %s\n',i,shortname);
        case 'Cancel'
            fprintf('Finished %i: %s\n',i,shortname);
            fprintf('... Stopping %s\n',mfilename);
            break;
    end
end

