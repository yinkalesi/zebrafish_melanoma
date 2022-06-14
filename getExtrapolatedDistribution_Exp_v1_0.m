%% getExtrapolatedDistribution_Exp_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 3/23/21
%  Project: Tumor Growth, Logarithmic Continuum Form

function [x_add,dist_add] = getExtrapolatedDistribution_Exp_v1_0(ntum,inc,size_range)

lx = log10(size_range);

d1_slope = -ntum/(lx(2)-lx(1));
ntum_add1 = round(ntum/inc);
delx = -inc/d1_slope;
x_add = 10.^((0:ntum_add1-1)*delx+lx(1))';
dist_add = ones(size(x_add))*inc;

end