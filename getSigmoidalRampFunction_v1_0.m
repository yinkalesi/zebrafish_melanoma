function [ramper] = getSigmoidalRampFunction_v1_0(center,width)
%% getSigmoidalRampFunction_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/12/17
%  Project: Tumor Growth, Logarithmic Continuum Form

%  get a function that evaluates to a sigmoidal curve that goes from 0.01 
%  to 0.99 from  x=center-width/2 to x=center+width/2

% (1+tanh(2.29756))/2 = 0.99
if(width>0)
    ramper = @(x) 0.5*(1+tanh((x-center)*2.29756*2/width));
else
    ramper = @(x) double(x>=center);
end