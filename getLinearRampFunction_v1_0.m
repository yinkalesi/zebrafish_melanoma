function [ramper] = getLinearRampFunction_v1_0(center,width)
%% getLinearRampFunction
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/12/17
%  Project: Tumor Growth, Logarithmic Continuum Form

%  get a function that evaluates to a ramp that goes from zero to one from
%  x=center-width/2 to x=center+width/2

if(width>0)
    ramper = @(x) max(0,min(1,0.5+(x-center)/width));
else
    ramper = @(x) double(x>=center);
end