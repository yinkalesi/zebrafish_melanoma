function [obj] = optimizer_v10_0_nograd(x0,key,fitter,chain)
%% optimizer_v10_0_nograd
%  Version 10.0
%  Author: Adeyinka Lesi
%  Date: 5/19/20

%  function obj = optimizer_v5_0(x0,key,data)
%  x0: 1x4 mattrix with log10 of growth and metastasis parameters
%  key: struct of parameters
%  (no longer used) data: struct for experimental data
%  fitter: struc holding fitting functions
%  chain: vector needed to calculate proper gradient
%% Version History
%  8.0: usng prePlotCalcs and gradient, and fit selection function for time
%  is max(fit)
%  8.1: changing fit selection criteria to log(e^fit/N) - this gives a
%  smooth function while maintaining the property that the worst fits are
%  heavily weighted
%  8.2: previously, forgot to factor in the fact that the log of x0 is
%  being used outside this function, which changes the proper gradient
%  value!
%  9.0: one more fitting parameter - carrying capacity
%  9.2: optimization for ramp parameters
%  9.3: conmod_v2_3 (changes to how initial tumor values are calculated)
%  9.4: using getFitter_v3_1 (with MMD fit func)
%  10.0: using conmod_3_1, getFitter_v4_0 & inferred data

% get parameters from pars
kg = x0(1:2);
kr = x0(3:4);
ks = x0([5 7]);
if(length(x0) < 8)
    % no carrying capacity, shed_meta_ratio, or reseeding parameters
    cc = inf;
    kc = [x0(6) inf 0 0];
    if(x0(6) ~= 0)
        kc(2) = x0(5)/x0(6);
    end
    ramp_center = 0;
    ramp_width = 0;
elseif(length(x0) < 11)
    % no shed_meta_ratio or reseeding parameters
    cc = x0(8);
    kc = [x0(6) inf 0 0];
    if(x0(6) ~= 0)
        kc(2) = x0(5)/x0(6);
    end
    ramp_center = x0(9);
    ramp_width = x0(10);
elseif(length(x0) < 12)
    % no reseeding parameters
    cc = x0(8);
    kc = [x0(6) x0(11) 0 0];
    ramp_center = x0(9);
    ramp_width = x0(10);
else
    cc = x0(8);
    kc = x0([6 11:13]);
    ramp_center = x0(9);
    ramp_width = x0(10);
end

changes = struct('GROWTH_PARAMETER1',kg(1),...
                 'GROWTH_PARAMETER2',kg(2),...
                 'DEATH_PARAMETER1',kr(1),...
                 'DEATH_PARAMETER2',kr(2),...
                 'SHED_PARAMETER1',ks(1),...
                 'SHED_PARAMETER2',ks(2),...
                 'CARRYING_CAPACITY',cc,...
                 'CTC_PARAMETER1',kc(1),...
                 'SHED_META_RATIO',kc(2),...
                 'CTC_PARAMETER2',(key.SHED_META_RATIO-1)*key.CTC_PARAMETER1,...
                 'SEED_PARAMETER1',kc(3),...
                 'SEED_PARAMETER2',kc(4),...
                 'DEATH_RAMP_CENTER',ramp_center,...
                 'DEATH_RAMP_WIDTH',ramp_width);
res = conmod_v3_1(key, changes); %conmod_v2_0 calculates cdfs as well

[fit0] = fitter.objective(res,res.key.inferred_data);
fit_weights = res.key.inferred_data.weights;
% fit_weights(1) = 0; % ignore fit of first time point
fit = fit0(fit_weights>0).*fit_weights(fit_weights>0); % ignore fit of first time point

% fitting factor (the lower the number, the less leeway for the fits at 
% different times to be different
q = 0.01;
exp_max = 300; % avoids returning inf
expfit = exp(min(fit/q,exp_max));
obj = sum(fit.*expfit)/sum(expfit);
if(isrow(chain))
    chain = chain';
end
% obj_grad = (grad*expfit'/sum(expfit)).*chain;

disp(['Parameters: ' res.key.TITLE]);
disp('Fit                              | Objective');
disp([num2str(fit0) ' | ' num2str(obj)]);
disp([num2str(fit_weights) ' | weights']);
% disp('Gradient | Objective');
% delim = ' | ';
% strline = delim;
% for i = 1:length(obj_grad)-1
%     strline = [strline;delim];
% end
% disp([num2str(grad) strline num2str(obj_grad)]);
disp(' ');

% Fit | Objective
% 0.028425      1.6703      1.0953 | 1.1349
% Gradient | Objective
% 0.001821214   -0.05190006    0.09281764 | 0.00038706
%    1.529625     -79.19107      141.9366 |    0.55284
end