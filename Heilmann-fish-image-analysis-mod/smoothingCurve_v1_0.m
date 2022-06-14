%% smoothingCurve_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 2/11/20
%% Version History
%  1.0: generate runing average of a curve

function [smooth] = smoothingCurve_v1_0(curve,nsmooth)
nc = length(curve);
nsh = floor(nsmooth/2);
smooth = curve;
navg = ones(size(curve));
for i = [1:nsh-1 nsh+1:nsmooth]
    starti = max(1,i-nsh+1);
    lasti = min(nc,i-nsh+1+nc);
    if(starti==1)
        smooth(starti:lasti) = smooth(starti:lasti)+curve(nc-lasti+1:end);
    else
        smooth(starti:lasti) = smooth(starti:lasti)+curve(1:lasti-starti+1);
    end
    navg(starti:lasti) = navg(starti:lasti)+1;
end
smooth = smooth./navg;

end