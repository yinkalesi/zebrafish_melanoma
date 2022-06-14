%% getTransformationRule_v1_4
%  Version 1.4
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 7/7/18
%  Project: Tumor Growth, Logarithmic Continuum Form
%  getTransformationRule returns a set of functions to guide the spacing of
%  the mesh in the simulation
%  key: struct; information on grid
%  transform: struct; contains transform, inverse transform, transform
%  derivative and double derivative

%% Version History
%  1.0: y=A*(x-1)+B*ln(x+C)
%  1.1: need to adjust for case where HARD_TUMOR_SIZE_LIMIT and
%  NUMBER_OF_SIZE_INTERVALS are close
%  1.2: adjustment to parameters so they work when xmax is moderately small
%  1.3: need to have decent guesses for zeros to find that
%  1.4: return null Transformation rule in certain cases (small xmax)
function [transform] = getTransformationRule_v1_4(key)

transform = struct();

xmax = key.HARD_TUMOR_SIZE_LIMIT;
Nc = key.NUMBER_OF_SIZE_INTERVALS;
r = min(100,xmax-1);
dy1 = 1;

zfunc = @(u) Nc.*(1-r/(xmax-1)).*log(abs(u)+1)-(1+Nc/(xmax-1)*(1-r/(xmax-1))).*log((xmax-1).*abs(u)+1);
guesses = 10.^(-10:0.1:0);
iguess = find(zfunc(guesses)>0,1);
if(isempty(iguess) || xmax-1 <= Nc) % xmax is small: no need for a transformation
    transform = getNullTransformation();
else
    tht = fzero(zfunc,guesses(iguess));
    if(isnan(tht))
        warning('Unable to find zero');
        tht = 0;
    else
        tht = abs(tht);
    end
    
    Ac = Nc*dy1*((xmax-1)*(r+Nc)-r*Nc)/(xmax-1)^3;
%     Ac = Nc*dy1/(xmax-1)*(1-(xmax-1-Nc)*(xmax-1-r)/(xmax-1)^2);
    Bc = (dy1-Ac)/log(tht+1);
    Cc = 1/tht-1;
    
    transform.x2y = @(x) Ac*(x-1)+Bc*log(x+Cc);
    if(Bc ~= 0)
        transform.y2x = @(y) lambertw(exp((y+Ac*(1+Cc))/Bc)*Ac/Bc)*Bc/Ac-Cc;
    else
        transform.y2x = @(y) y/Ac+1;
    end
    transform.dydx = @(x) Ac+Bc./(x+Cc);
    transform.d2ydx2 = @(x) -Bc./(x+Cc).^2;
    transform.params = struct('Ac',Ac,'Bc',Bc,'Cc',Cc);
end



