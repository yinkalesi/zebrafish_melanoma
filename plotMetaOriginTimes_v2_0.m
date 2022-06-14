%% plotMetaOriginTimes_v2_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Dr. Silja Heilmann
%  Date: 2/8/21
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  2.0: the last observed time for each metastasis is chosen as the point
% to backtrack from :<

function [] = plotMetaOriginTimes_v2_0(res)

meta_origin = zeros(size(res.key.META_ORIGIN(:,end)));
for j = 1:length(meta_origin)
    [i1,i2] = find(res.data.alternate(end).origin_time_loc==j,1,'first');
    if(res.data.alternate(end).origin_time_corr(i1,end)>=0)
        meta_origin(j) = res.data.alternate(end).origin_time_corr(i1,end);
    else
        meta_origin(j) = res.data.alternate(end).origin_time_corr(i1,i2);
    end
end
[~,mord]=sort(meta_origin);
orig_times = meta_origin(mord);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));

stepx = zeros(length(orig_times)*2,1);
stepcdf = zeros(length(orig_times)*2,1);

stepx((1:length(orig_times))*2-1) = orig_times;
stepx((1:length(orig_times))*2) = orig_times;
stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
stepcdf((1:length(orig_times))*2) = met_cumdist;

figure;
plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-','LineWidth',1.5)
legend('Model','Data');
xlabel('Time (Days)');
ylabel('Number of Metastasis (Scaled)');
end