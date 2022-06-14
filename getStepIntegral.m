%% getStepIntegral
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 8/6/15
%  Project: Tumor Growth, Logarithmic Continuum Form

function integral = getStepIntegral(x,y)
% getStepIntegral adds up the ys at and before each x value and stores that
% in the integral vector
% x = 1xM vector
% y = 1xM vector
% integral = 1xM vector

integral = zeros(size(x));
integral(1) = y(1);
for i = 2:length(x)
    integral(i) = integral(i-1)+y(i);
end