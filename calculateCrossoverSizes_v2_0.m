%% calculateCrossoverSizes_v2_0.m
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 1/11/22
%% Version History
%  1.1: for MM parameters
%  1.2: made into function; uses key as input and assumes either MM11,
%  MM123 or PL form for reduction rates

function [xbals] =  calculateCrossoverSizes_v2_0(key)

g1 = key.PARAMETERS.growth1;
g2 = key.PARAMETERS.growth2;
d1 = key.PARAMETERS.death1;
d2 = key.PARAMETERS.death2;
s1 = key.PARAMETERS.shed1;
s2 = key.PARAMETERS.shed2;

if(isequal(key.RATE_GENERATOR,@getRates_michaelisMenten_form123_v1_0))
    xb0 = (d1/d2/g1)^(3/(3*g2-1)); % guess (can use solver to get closer)
    pfun = @(x) (1+d2*x^(2/3))*(1-s1/g1*x^(s2-g2))-d1/g1*x^(1-g2);
    
elseif(isequal(key.RATE_GENERATOR,@getRates_michaelisMenten_form11_v1_0))
    xb0 = (d1/d2/g1)^(1/g2);
    pfun = @(x) (1+d2*x)*(1-s1/g1*x^(s2-g2))-d1/g1*x^(1-g2);
else
   % assuming power law 
   xb0 = (d1/g1)^(1/(g2-d2));
   pfun = @(x) 1-d1/g1*x^(d2-g2)-s1/g1*x^(s2-g2);
end

xbals = fzero(pfun,xb0);
