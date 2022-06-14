%% getSpatialThreshold
%  Version 1.2
%  Author: Adeyinka Lesi
%  Date: 2/12/20
%  Project: Fish Image Analysis
%  find the best position (index) to separate an array of data so the two
%  sides have maximally different inter-class variance
%% Version History
%  1.1: (2/12/20) A reminder that the (cands-1).*(length(data)-cands+1)
%  weighting term was added to favor leaving the threshold near the center.
%  This is because this method is used for cases that aren't clean cut; in
%  those case the obj function peak would be near an edge without applying
%  the weight. The literature backing for this technique is that this is
%  basically a version of Otsu's method
%  1.2: to fix issue where there are two potential thresholds are available
%  (obj has two peaks). Will force data to be monotonically increasing or
%  decreasing which has the effect of a) smoothing out the data hopefully
%  and b) eliminating any dips in the data that might lead to errors

function thresh = getSpatialThreshold_v1_2(data_in)

data = data_in;
end_depth = ceil(length(data)/10);
left_end = mean(data(1:end_depth));
right_end = mean(data(end-end_depth+1:end));

data(1) = left_end;
data(end) = right_end;
if(left_end<right_end)
    % increasing
    for i = 2:length(data)
        if(data(i)<data(i-1))
            data(i)=data(i-1);
        end
    end    
else
    % decreasing
    for i = 2:length(data)
        if(data(i)>data(i-1))
            data(i)=data(i-1);
        end
    end
end

cands = 2:length(data);
left_avg = zeros(size(cands));
right_avg = zeros(size(cands));

for i = 1:length(cands)
    left_avg(i) = mean(data(1:cands(i)-1));
    right_avg(i) = mean(data(cands(i):end));
end

obj = (cands-1).*(length(data)-cands+1).*(left_avg-right_avg).^2;
[~,imax] = max(obj);
thresh = cands(imax);

% figure;
% plot(1:length(data),data/sum(data),1:length(data_in),data_in/sum(data_in),'--',cands,obj/sum(obj),cands,(left_avg-right_avg).^2/sum((left_avg-right_avg).^2),'--');
end