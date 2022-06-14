%% plotMetaOriginTimes_v2_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Dr. Silja Heilmann
%  Date: 2/12/21
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  2.0: the last observed time for each metastasis is chosen as the point
% to backtrack from :<
%  2.1: plotting both forms of metastasis tracking

function [] = plotMetaOriginTimes_v2_1(res)

meta_origin1 = res.key.META_ORIGIN(:,end);
meta_origin2 = zeros(size(res.key.META_ORIGIN(:,end)));
for j = 1:length(meta_origin2)
    [i1,i2] = find(res.data.alternate(end).origin_time_loc==j,1,'first');
    if(res.data.alternate(end).origin_time_corr(i1,end)>=0)
        meta_origin2(j) = res.data.alternate(end).origin_time_corr(i1,end);
    else
        meta_origin2(j) = res.data.alternate(end).origin_time_corr(i1,i2);
    end
end
[~,mord1]=sort(meta_origin1);
orig_times1 = meta_origin1(mord1);
met_cumdist1 = cumsum(res.key.META_ORIGIN(mord1,end-1));
[~,mord2]=sort(meta_origin2);
orig_times2 = meta_origin2(mord2);
met_cumdist2 = cumsum(res.key.META_ORIGIN(mord2,end-1));

% calculate third compromise tumor count
min1 = min(orig_times1);
lmind1 = find(orig_times1==min1,1,'last');
min2 = min(orig_times2);
lmind2 = find(orig_times2==min2,1,'last');

if(min2>min1 || (min2==min1&&lmind2<lmind1))
    % append vestigial tumors from orig_times1 as extra tumors to
    % orig_times2 to make orig_times3
    inter = find(orig_times1==min2);
    if(isempty(inter))
        junc = find(orig_times1<min2,1,'last');
    else
        shift = max(0,length(inter)-lmind2);
        junc = inter(1)+shift-1;
    end
    orig_times3 = [orig_times1(1:junc); orig_times2];
    met_cumdist3 = [met_cumdist1(1:junc); met_cumdist1(junc)+met_cumdist2];
elseif(min2<min1 || (min2==min1&&lmind2>lmind1))
    % append vestigial tumors from orig_times2 as extra tumors to
    % orig_times1 to make orig_times3
    inter = find(orig_times2==min1);
    if(isempty(inter))
        junc = find(orig_times2<min1,1,'last');
    else
        shift = max(0,length(inter)-lmind1);
        junc = inter(1)+shift-1;
    end
    orig_times3 = [orig_times2(1:junc); orig_times1];
    met_cumdist3 = [met_cumdist2(1:junc); met_cumdist2(junc)+met_cumdist1];
else
    orig_times3 = orig_times2;
    met_cumdist3 = met_cumdist2;
end

stepx1 = zeros(length(orig_times1)*2,1);
stepcdf1 = zeros(length(orig_times1)*2,1);

stepx1((1:length(orig_times1))*2-1) = orig_times1;
stepx1((1:length(orig_times1))*2) = orig_times1;
stepcdf1((1:length(orig_times1))*2-1) = [0; met_cumdist1(1:end-1)];
stepcdf1((1:length(orig_times1))*2) = met_cumdist1;

stepx2 = zeros(length(orig_times2)*2,1);
stepcdf2 = zeros(length(orig_times2)*2,1);

stepx2((1:length(orig_times2))*2-1) = orig_times2;
stepx2((1:length(orig_times2))*2) = orig_times2;
stepcdf2((1:length(orig_times2))*2-1) = [0; met_cumdist2(1:end-1)];
stepcdf2((1:length(orig_times2))*2) = met_cumdist2;

stepx3 = zeros(length(orig_times3)*2,1);
stepcdf3 = zeros(length(orig_times3)*2,1);

stepx3((1:length(orig_times3))*2-1) = orig_times3;
stepx3((1:length(orig_times3))*2) = orig_times3;
stepcdf3((1:length(orig_times3))*2-1) = [0; met_cumdist3(1:end-1)];
stepcdf3((1:length(orig_times3))*2) = met_cumdist3;

figure;
plot(res.t3,cumsum(res.num_meta),stepx1,stepcdf1,'o-',stepx2,stepcdf2,'s-',stepx3,stepcdf3,'d-','LineWidth',1.5)
legend('Model','Data (Detect Time)','Data (Last Time)','Data (Combo Time)');
xlabel('Time (Days)');
ylabel('Number of Metastasis (Scaled)');
end