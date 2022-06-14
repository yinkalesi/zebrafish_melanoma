%% printParameters_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 7/29/20
%  Project: Tumor Growth, Logarithmic Continuum Form

function [par_str,format] = printParameters_v1_0(pars,key)

if(pars(9) > 0)
    ramp_start = pars(9)-0.5*pars(10);
    ramp_end = pars(9)+0.5*pars(10);
    ramp_start = round(10*ramp_start)*0.1;
    ramp_end = round(10*ramp_end)*0.1;
    format = ['k_g=%1.2e x^{%1.2f}, ' ...
        'k_r=%1.2e x^{%1.2f} (ramp(' ...
        num2str(ramp_start) '->' num2str(ramp_end) ')), '...
        '(k_s=%1.2e)k_m=%1.2e x^{%1.2f}, ' ...
        'cc=%1.2e'];
else
    format = ['k_g=%1.2e x^{%1.2f}, ' ...
        'k_r=%1.2e x^{%1.2f}, ' ...
        '(k_s=%1.2e)k_m=%1.2e x^{%1.2f}, ' ...
        'cc=%1.2e'];
end

par_str = sprintf(format,pars(1:8));

% fprintf('%s\n',par_str);
