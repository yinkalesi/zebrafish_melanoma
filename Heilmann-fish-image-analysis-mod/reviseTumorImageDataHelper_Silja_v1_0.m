%% reviseTumorImageDataHelper_Silja_v1_0
%  version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 3/11/21
%% Version History
%  1.0: from reviseTumorImageDataHelper_v1_2, for Silja's data

dat_dir = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\Continuous Tumor Growth Model\Heilmann-fish-image-analysis-mod\'; 
% batch1 = {...
%     [basedir 'Silja_022520_r2_Batch1_DORSAL_1e5_summary'],...
%     [basedir 'Silja_022520_r2_Batch1_DORSAL_1e6_summary'],...
%     [basedir 'Silja_022520_r2_Batch1_DORSAL_5e5_summary'],...
%     [basedir 'Silja_022520_r2_Batch1_VENTRAL_1e5_summary'],...
%     [basedir 'Silja_022520_r2_Batch1_VENTRAL_1e6_summary'],...
%     [basedir 'Silja_022520_r2_Batch1_VENTRAL_5e5_summary']};
% 
% batch2 = {...
%     [basedir 'Silja_022520_r2_Batch2_DORSAL_1e5_summary'],...
%     [basedir 'Silja_022520_r2_Batch2_DORSAL_1e6_summary'],...
%     [basedir 'Silja_022520_r2_Batch2_DORSAL_5e5_summary'],...
%     [basedir 'Silja_022520_r2_Batch2_VENTRAL_1e5_summary'],...
%     [basedir 'Silja_022520_r2_Batch2_VENTRAL_1e6_summary'],...
%     [basedir 'Silja_022520_r2_Batch2_VENTRAL_5e5_summary']};

summary_files = dir([dat_dir 'Silja_022520_r2*tumors.mat']);
image_files =  dir([dat_dir 'Silja_022520_r2*images.mat']);

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
    if(~exist([newname '_tumors_revised_test.mat'],'file'))
        reviseTumorImageData_v1_2(tum_file,img_file,newname);
    else
        fprintf('File %s already exists\n',[newname '_tumors_revised_test.mat']);
    end
    move_on = questdlg(['Move on with ' shortname '?']);
    switch move_on
        case 'Yes'
            fprintf('Finished %i: %s\n',i,shortname);
        case 'No'
            fprintf('Redoing %i: %s\n',i,shortname);
            reviseTumorImageData_v1_2(tum_file,img_file,newname,[newname '_tumors_revised_test.mat']);
            fprintf('Finished %i: %s\n',i,shortname);
        case 'Cancel'
            fprintf('Finished %i: %s\n',i,shortname);
            fprintf('... Stopping %s\n',mfilename);
            break;
    end
end

