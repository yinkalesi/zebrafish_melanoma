%% getRates_michaelisMenten_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 11/10/20
%  Project: Tumor Growth, Logarithmic Continuum Form
% getRates_michaelisMenten_v1_0 uses key.PARAMETERS to create functions for
% the growth, death, shedding and metastasis of tumors. 
% params: struct; rate parameters
% rate_funcs: struct; contains rate functions
% string_form: string; the formula for parameters is stated explicitly
%% VERSION HISTORY
%  1.0: from getRates_powerLaw_v5_0, This version assume
% growth is described as k*j^i and death as form k1*x/(1+k2*x)

function [rate_funcs,string_form] = getRates_michaelisMenten_form11_v1_0(key)

params = key.PARAMETERS;
rate_funcs = struct();
rate_funcs.growth = @(x) params.growth1*x.^params.growth2;
rate_funcs.shed = @(x) params.shed1*x.^params.shed2.*(x>=2);
rate_funcs.seed = @(x) params.seed1*x.^params.seed2;
rate_funcs.meta_ss = @(x,sf) params.ctc1*rate_funcs.shed(x)./(sf.seed+params.ctc1+params.ctc2);
rate_funcs.meta = @(x) rate_funcs.shed(x)./(params.shed_meta_ratio); % for backward compatibility; try not to use this
uncut_func = @(x) params.death1*x./(1+params.death2*x);
if(key.USING_DEATH_CUTOFF)
    rate_funcs.death = @(x) uncut_func(x).*(x<params.death_cutoff);
else
    rate_funcs.death = uncut_func;
end
% first derivative
rate_funcs.growth_deriv = @(x) params.growth2*params.growth1*x.^(params.growth2-1);
rate_funcs.shed_deriv = @(x) params.shed2*params.shed1*x.^(params.shed2-1).*(x>=2);
rate_funcs.seed_deriv = @(x) params.seed2*params.seed1*x.^(params.seed2-1);
uncut_func2 = @(x) params.death1./(1+params.death2*x).^2;
dcv = 0.1;
gfunc = @(x,x0) -exp(-0.5*(x-x0).^2/dcv)/dcv/sqrt(2*pi);
gfunc_deriv = @(x,x0) (x-x0).*exp(-0.5*(x-x0).^2/dcv)/dcv^3/sqrt(2*pi);
if(key.USING_DEATH_CUTOFF)
    rate_funcs.death_deriv = @(x) uncut_func2(x).*(x<params.death_cutoff)...
        +uncut_func(x).*gfunc(x,params.death_cutoff);
else
    rate_funcs.death_deriv = uncut_func2;
end
% second derivative
rate_funcs.growth_2deriv = @(x) params.growth2*(params.growth2-1)*params.growth1*x.^(params.growth2-2);
rate_funcs.shed_2deriv = @(x) params.shed2*(params.shed2-1)*params.shed1*x.^(params.shed2-2).*(x>=2);
rate_funcs.seed_2deriv = @(x) params.seed2*(params.seed2-1)*params.seed1*x.^(params.seed2-2);
uncut_func3 = @(x) -2*params.death2*params.death1./(1+params.death2*x).^3;
if(key.USING_DEATH_CUTOFF)
    rate_funcs.death_2deriv = @(x) uncut_func3(x).*(x<params.death_cutoff)...
        +2*uncut_func2(x).*gfunc(x,params.death_cutoff)...
        +uncut_func(x).*gfunc_deriv(x,params.death_cutoff);
else
    rate_funcs.death_2deriv = uncut_func3;
end

% parameter derivatives
rate_funcs.growth_parderiv = @(x,parname) getGrowthParDeriv(x,parname,params,key);
rate_funcs.death_parderiv = @(x,parname) getDeathParDeriv(x,parname,params,key);
rate_funcs.shed_parderiv = @(x,parname) getShedParDeriv(x,parname,params,key);
rate_funcs.meta_parderiv = @(x,parname) getShedParDeriv(x,parname,params,key)./(params.shed_meta_ratio); % for backward compatibility; try not to use this
rate_funcs.seed_parderiv = @(x,parname) getSeedParDeriv(x,parname,params,key);
rate_funcs.size1_ratio = params.growth1/(params.growth1+params.death1/(1+params.death2));
rate_funcs.size1_ratio_parderiv = @(parname) getSize1RatioDeriv(parname,params,key);
rate_funcs.meta_ss_parderiv = @(x,parname,sf) getMetaSSParDeriv(x,parname,params,key,sf);

% string containing the form of the equations
if(key.USING_DEATH_CUTOFF)
    if(isfield(key,'USING_DEATH_STARTER') && key.USING_DEATH_STARTER)
        format = ['k_g=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_r=%1.2e\\cdot{}x/(1+%1.2e\\cdot{}x) (x<' ...
            num2str(params.death_cutoff) ', t>' ...
            num2str(params.death_start_time) '), ' ...
            'k_s=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_0=%1.2e, k_1=%1.2e, ' ...
            'k_2=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'cc=%1.2e'];
    elseif(isfield(key,'USING_DEATH_RAMP') && key.USING_DEATH_RAMP)
        ramp_start = params.death_ramp_center-0.5*params.death_ramp_width;
        ramp_end = params.death_ramp_center+0.5*params.death_ramp_width;
        ramp_start = round(10*ramp_start)*0.1;
        ramp_end = round(10*ramp_end)*0.1;
        format = ['k_g=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_r=%1.2e\\cdot{}x/(1+%1.2e\\cdot{}x) (x<' ...
            num2str(params.death_cutoff) ',ramp(' ...
            num2str(ramp_start) '->' num2str(ramp_end) ')), '...
            'k_s=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_0=%1.2e, k_1=%1.2e, ' ...
            'k_2=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'cc=%1.2e'];
    else
        format = ['k_g=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_r=%1.2e\\cdot{}x/(1+%1.2e\\cdot{}x) (x<' ...
            num2str(params.death_cutoff) '), ' ...
            'k_s=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_0=%1.2e, k_1=%1.2e, ' ...
            'k_2=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'cc=%1.2e'];
    end
else
    if(isfield(key,'USING_DEATH_STARTER') && key.USING_DEATH_STARTER)
        format = ['k_g=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_r=%1.2e\\cdot{}x/(1+%1.2e\\cdot{}x) (t>' ...
            num2str(params.death_start_time) '), ' ...
            'k_s=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_0=%1.2e, k_1=%1.2e, ' ...
            'k_2=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'cc=%1.2e'];
    elseif(isfield(key,'USING_DEATH_RAMP') && key.USING_DEATH_RAMP)
        ramp_start = params.death_ramp_center-0.5*params.death_ramp_width;
        ramp_end = params.death_ramp_center+0.5*params.death_ramp_width;
        ramp_start = round(10*ramp_start)*0.1;
        ramp_end = round(10*ramp_end)*0.1;
        format = ['k_g=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_r=%1.2e\\cdot{}x/(1+%1.2e\\cdot{}x) (ramp(' ...
            num2str(ramp_start) '->' num2str(ramp_end) ')), '...
            'k_s=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_0=%1.2e, k_1=%1.2e, ' ...
            'k_2=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'cc=%1.2e'];
    else
        format = ['k_g=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_r=%1.2e\\cdot{}x/(1+%1.2e\\cdot{}x), ' ...
            'k_s=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'k_0=%1.2e, k_1=%1.2e, ' ...
            'k_2=%1.2e\\cdot{}x^{%1.2f}, ' ...
            'cc=%1.2e'];
    end
end

string_form = sprintf(format,[params.growth1, params.growth2, ...
                              params.death1, params.death2, ...
                              params.shed1, params.shed2, ...
                              params.ctc1, params.ctc2, ...
                              params.seed1, params.seed2, ...
                              params.carrying_capacity]);

end

function [pardevs,iszero] = getGrowthParDeriv(x,parname,pars,key)
iszero = true;
switch parname
    case 'growth1'
        pardevs = x.^(pars.growth2);
        iszero = false;
    case 'growth2'
        pardevs = pars.growth1*x.^(pars.growth2).*log(x);
        iszero = false;
    case 'death1'
        pardevs = zeros(size(x));
    case 'death2'
        pardevs = zeros(size(x));
    case 'shed1'
        pardevs = zeros(size(x));
    case 'shed2'
        pardevs = zeros(size(x));
    case 'seed1'
        pardevs = zeros(size(x));
    case 'seed2'
        pardevs = zeros(size(x));
    case 'ctc1'
        pardevs = zeros(size(x));
    case 'ctc2'
        pardevs = zeros(size(x));
    case 'carrying_capacity'
        pardevs = zeros(size(x));
    case 'death_ramp_center'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_ramp_width'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_cutoff'
        if(key.USING_DEATH_CUTOFF)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_start_time'
        if(key.USING_DEATH_STARTER)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    otherwise
        pardevs = zeros(size(x));
end
end

function [pardevs,iszero] = getDeathParDeriv(x,parname,pars,key)
iszero = true;
switch parname
    case 'growth1'
        pardevs = zeros(size(x));
    case 'growth2'
        pardevs = zeros(size(x));
    case 'death1'
        pardevs = x./(1+pars.death2*x);
        iszero = false;
    case 'death2'
        pardevs = -pars.death1*x.^2./(1+pars.death2*x).^2;
        iszero = false;
    case 'shed1'
        pardevs = zeros(size(x));
    case 'shed2'
        pardevs = zeros(size(x));
    case 'seed1'
        pardevs = zeros(size(x));
    case 'seed2'
        pardevs = zeros(size(x));
    case 'ctc1'
        pardevs = zeros(size(x));
    case 'ctc2'
        pardevs = zeros(size(x));
    case 'carrying_capacity'
        pardevs = zeros(size(x));
    case 'death_ramp_center'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_ramp_width'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_cutoff'
        if(key.USING_DEATH_CUTOFF)
            % derivative is nonzero close to cutoff size but zero
            % everywhere else; will use a gaussian distribution to
            % approximate this
            mdif = min(abs(x-pars.death_cutoff));
            Ddc = max(mdif*0.5,0.5);
            pardevs = exp(-(x-pars.death_cutoff).^2/(4*Ddc))/sqrt(pi*Ddc)*0.5;
            iszero = false;
        else
            pardevs = zeros(size(x));
        end
    case 'death_start_time'
        if(key.USING_DEATH_STARTER)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    otherwise
        pardevs = zeros(size(x));
end
end

function [pardevs,iszero] = getShedParDeriv(x,parname,pars,key)
iszero = true;
switch parname
    case 'growth1'
        pardevs = zeros(size(x));
    case 'growth2'
        pardevs = zeros(size(x));
    case 'death1'
        pardevs = zeros(size(x));
    case 'death2'
        pardevs = zeros(size(x));
    case 'shed1'
        pardevs = x.^(pars.shed2);
        iszero = false;
    case 'shed2'
        pardevs = pars.shed1*x.^(pars.shed2).*log(x);
        iszero = false;
    case 'seed1'
        pardevs = zeros(size(x));
    case 'seed2'
        pardevs = zeros(size(x));
    case 'ctc1'
        pardevs = zeros(size(x));
    case 'ctc2'
        pardevs = zeros(size(x));
    case 'carrying_capacity'
        pardevs = zeros(size(x));
    case 'death_ramp_center'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_ramp_width'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_cutoff'
        if(key.USING_DEATH_CUTOFF)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_start_time'
        if(key.USING_DEATH_STARTER)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    otherwise
        pardevs = zeros(size(x));
end
end

function [pardevs,iszero] = getSeedParDeriv(x,parname,pars,key)
iszero = true;
switch parname
    case 'growth1'
        pardevs = zeros(size(x));
    case 'growth2'
        pardevs = zeros(size(x));
    case 'death1'
        pardevs = zeros(size(x));
    case 'death2'
        pardevs = zeros(size(x));
    case 'shed1'
        pardevs = zeros(size(x));
    case 'shed2'
        pardevs = zeros(size(x));
    case 'seed1'
        pardevs = x.^(pars.seed2);
        iszero = false;
    case 'seed2'
        pardevs = pars.seed1*x.^(pars.seed2).*log(x);
        iszero = false;
    case 'ctc1'
        pardevs = zeros(size(x));
    case 'ctc2'
        pardevs = zeros(size(x));
    case 'carrying_capacity'
        pardevs = zeros(size(x));
    case 'death_ramp_center'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_ramp_width'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_cutoff'
        if(key.USING_DEATH_CUTOFF)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_start_time'
        if(key.USING_DEATH_STARTER)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    otherwise
        pardevs = zeros(size(x));
end
end

function [pardevs,iszero] = getSize1RatioDeriv(parname,pars,key)
iszero = true;
switch parname
    case 'growth1'
        pardevs = pars.death1/(pars.growth1+pars.death1/(1+pars.death2))^2;
        iszero = false;
    case 'growth2'
        pardevs = 0;
    case 'death1'
        pardevs = -pars.growth1/(1+pars.death2)/(pars.growth1+pars.death1/(1+pars.death2))^2;
        iszero = false;
    case 'death2'
        pardevs = pars.death1*pars.growth1/(pars.growth1*(1+pars.death2)+pars.death1)^2;
    case 'shed1'
        pardevs = 0;
    case 'shed2'
        pardevs = 0;
    case 'seed1'
        pardevs = 0;
    case 'seed2'
        pardevs = 0;
    case 'ctc1'
        pardevs = 0;
    case 'ctc2'
        pardevs = 0;
    case 'carrying_capacity'
        pardevs = 0;
    case 'death_ramp_center'
        if(key.USING_DEATH_RAMP)
            pardevs = 0;
        else
            pardevs = 0;
        end
    case 'death_ramp_width'
        if(key.USING_DEATH_RAMP)
            pardevs = 0;
        else
            pardevs = 0;
        end
    case 'death_cutoff'
        if(key.USING_DEATH_CUTOFF)
            pardevs = 0;
        else
            pardevs = 0;
        end
    case 'death_start_time'
        if(key.USING_DEATH_STARTER)
            pardevs = 0;
        else
            pardevs = 0;
        end
    otherwise
        pardevs = 0;
end
end

function [pardevs,iszero] = getMetaSSParDeriv(x,parname,pars,key,sf)
iszero = true;
switch parname
    case 'growth1'
        pardevs = zeros(size(x));
    case 'growth2'
        pardevs = zeros(size(x));
    case 'death1'
        pardevs = zeros(size(x));
    case 'death2'
        pardevs = zeros(size(x));
    case 'shed1'
        pardevs = pars.ctc1./(sf.seed+pars.ctc1+pars.ctc2).*x.^(pars.shed2);
        iszero = false;
    case 'shed2'
        pardevs = pars.ctc1./(sf.seed+pars.ctc1+pars.ctc2)*pars.shed1.*x.^(pars.shed2).*log(x);
        iszero = false;
    case 'seed1'
        pardevs = -pars.ctc1./(sf.seed+pars.ctc1+pars.ctc2).^2.*sf.seed_deriv1*pars.shed1.*x.^(pars.shed2);
        iszero = false;
    case 'seed2'
        pardevs = -pars.ctc1./(sf.seed+pars.ctc1+pars.ctc2).^2.*sf.seed_deriv2*pars.shed1.*x.^(pars.shed2);
        iszero = false;
    case 'ctc1'
        pardevs = (1./(sf.seed+pars.ctc1+pars.ctc2)-pars.ctc1./(sf.seed+pars.ctc1+pars.ctc2).^2)*pars.shed1.*x.^(pars.shed2);
        iszero = false;
    case 'ctc2'
        pardevs = -pars.ctc1./(sf.seed+pars.ctc1+pars.ctc2).^2*pars.shed1.*x.^(pars.shed2);
        iszero = false;
    case 'carrying_capacity'
        pardevs = zeros(size(x));
    case 'death_ramp_center'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_ramp_width'
        if(key.USING_DEATH_RAMP)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_cutoff'
        if(key.USING_DEATH_CUTOFF)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    case 'death_start_time'
        if(key.USING_DEATH_STARTER)
            pardevs = zeros(size(x));
        else
            pardevs = zeros(size(x));
        end
    otherwise
        pardevs = zeros(size(x));
end
end
