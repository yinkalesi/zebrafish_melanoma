%% plotMetaOriginTimes_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Dr. Silja Heilmann
%  Date: 6/2/20
%  Project: Tumor Growth, Logarithmic Continuum Form

function [] = plotMetaOriginTimes_v1_0(res)

[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
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