function [clean] = removeOverlappingRegions_v1_0(bw,flaws,thresh)
%% removeOverlappingRegions_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 9/25/17
%  Project: Fish Image Analysis
%  remove known bad signal completely from a bw image (simply subtracting)
%  doesn't account for the sizes of the regions in the bw and the flaws
%  images

regs = bwlabel(bw);
nregs = max(regs(:));
clean = false(size(bw));

% check each region for flaws
for r = 1:nregs
    flaw_count = sum(flaws(regs==r));
    reg_size = sum(bw(regs==r));
    
    if(flaw_count/reg_size < thresh)
        clean(regs==r) = 1;
    end
end