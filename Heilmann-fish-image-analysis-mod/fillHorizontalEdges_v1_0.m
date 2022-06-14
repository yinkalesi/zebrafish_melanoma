function [bw] = fillHorizontalEdges_v1_0(xneg,xpos)
% xneg is non-positive and xpos is non-negative
    bw = zeros(size(xneg));
    for row = 1:size(bw,1)
        bw(row,:) = markValleyToPeak(xneg(row,:),xpos(row,:));
    end
end

function [marked] = markValleyToPeak(valleys,peaks)
% input is a non-positive curve (valleys) and a non-negative curve (peaks)
% The idea is to find the extrema (vallyes and peaks) of section of curve 
% (a section is surrounded be zeros on either end). marked will be a binary
% line where regions with a valley to the left and a peak to the right are
% marked with ones and everything else is a zero

% 1) find valleys and peaks
valley_starts = find((valleys(1:end-1)==0).*((valleys(2:end)-valleys(1:end-1))~=0));
valley_locs = zeros(1,length(valley_starts));
for i = 1:length(valley_starts)-1
    [~,loc] = min(valleys(valley_starts(i):valley_starts(i+1)-1));
    valley_locs(i) = loc-1+valley_starts(i);
end
% do final iteration separately since range is to end of line
vlen = length(valley_starts);
if(vlen > 0)
    [~,loc] = min(valleys(valley_starts(vlen):end));
    valley_locs(vlen) = loc-1+valley_starts(vlen);
end

peak_starts = find((peaks(1:end-1)==0).*((peaks(2:end)-peaks(1:end-1))~=0));
peak_locs = zeros(1,length(peak_starts));
for i = 1:length(peak_starts)-1
    [~,loc] = max(peaks(peak_starts(i):peak_starts(i+1)-1));
    peak_locs(i) = loc-1+peak_starts(i);
end
% do final iteration separately since range is to end of line
plen = length(peak_starts);
if(plen > 0)
    [~,loc] = max(peaks(peak_starts(plen):end));
    peak_locs(plen) = loc-1+peak_starts(plen);
end

% 2) pair up valleys and peaks to make marked
marked = zeros(1,length(valleys));
if(vlen > 0)
    next_valley = 1;
    next_peak = find(peak_locs>valley_locs(next_valley),1);
    while(next_valley < vlen && ~isempty(next_peak))
        if(valley_locs(next_valley+1) >= peak_locs(next_peak))
            % case 1: the closest extrema after the current valley is a peak
            marked(valley_locs(next_valley):peak_locs(next_peak)) = 1;
        end
        % case 2: the closest extrema after the current valley is a valley and no marking needed
        % iterate
        next_valley = next_valley+1;
        next_peak = find(peak_locs>valley_locs(next_valley),1);
    end
    % handle final valley as long as vlen > 0 and there is another peak
    if(~isempty(next_peak))
        % variables have already been iterated
        marked(valley_locs(next_valley):peak_locs(next_peak)) = 1;
    end
end
    
end

