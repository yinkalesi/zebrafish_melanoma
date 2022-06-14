%% plotMetaOriginTimes_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Dr. Silja Heilmann
%  Date: 6/2/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  1.1: plot two metastasis plots

function [] = plotMetaOriginTimes_v1_1(res1,res2,desc1,desc2)

[~,mord1]=sort(res1.key.META_ORIGIN(:,end));
orig_times1 = res1.key.META_ORIGIN(mord1,end);
met_cumdist1 = cumsum(res1.key.META_ORIGIN(mord1,end-1));

stepx1 = zeros(length(orig_times1)*2,1);
stepcdf1 = zeros(length(orig_times1)*2,1);

stepx1((1:length(orig_times1))*2-1) = orig_times1;
stepx1((1:length(orig_times1))*2) = orig_times1;
stepcdf1((1:length(orig_times1))*2-1) = [0; met_cumdist1(1:end-1)];
stepcdf1((1:length(orig_times1))*2) = met_cumdist1;

[~,mord2]=sort(res2.key.META_ORIGIN(:,end));
orig_times2 = res2.key.META_ORIGIN(mord2,end);
met_cumdist2 = cumsum(res2.key.META_ORIGIN(mord2,end-1));

stepx2 = zeros(length(orig_times2)*2,1);
stepcdf2 = zeros(length(orig_times2)*2,1);

stepx2((1:length(orig_times2))*2-1) = orig_times2;
stepx2((1:length(orig_times2))*2) = orig_times2;
stepcdf2((1:length(orig_times2))*2-1) = [0; met_cumdist2(1:end-1)];
stepcdf2((1:length(orig_times2))*2) = met_cumdist2;

figure;
plot(res1.t3,cumsum(res1.num_meta),'b-',stepx1,stepcdf1,'bo--',...
    res2.t3,cumsum(res2.num_meta),'r-',stepx2,stepcdf2,'rs--','LineWidth',1.5)
legend({[desc1 ', Model'],[desc1 ', Data'],[desc2 ', Model'],[desc2 ', Data']},'Location','Best');
xlabel('Time (Days)');
ylabel('Number of Metastasis (Scaled)');
end