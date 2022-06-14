%% runModelFunction_v1_1
%  Version 1.1
%  Author: Adeyinka Lesi
%  Date: 3/26/18
%  Project: Tumor Growth, Logarithmic Continuum Form

function [mg_res,data,mg_fit] = runModelFunction_v1_1(datafile,params,conv,select)

% set name
% name by combinining data file with time
slashloc = find(datafile=='/',1,'last');
if(isempty(slashloc))
    slashloc = 0;
end
dotloc = find(datafile=='.',1,'last');
if(isempty(dotloc))
    dotloc = length(datafile)+1;
end
clockdat = clock;
name = [datafile(slashloc+1:dotloc-1) '_' sprintf('%02i%02i%4i-%02i%02i',clockdat([2 3 1 4 5]))];

if(~exist('conv','var'))
    conv = getNullConverter();
end
if(~exist('select','var'))
    select = [];
end

% load data
if(~isempty(select))
    data = getExperimentalData_v3_0(datafile,conv,select);
else
    data = getExperimentalData_v3_0(datafile,conv);
end
data.plot_indices = 1:6;
data.cdf_adj = [0 0 0 0 1 1];

% initial guess
kg = params(1:2);
kr = params(3:4);
ks = params(5);
km = params(6:7);
cc = params(8);
ramp_center = params(9);
ramp_width = params(10);
use_ramp = ramp_center > 0;

% key = getKey_metagen_v1_1([79.616798681379152   0.569437688675373],[0.182851433668764   0.947026744394068],ks,[ks*1e-3 gamma], ...
%     data.t(end),data.t(end)*10,data.x,data.dist(:,1),2000,primary_init,data.lowx,name);
params = [kg kr ks km];
key = getDataKey_v3_0(data,1000,1,name,params,'uniform');

changes = struct('SOLUTION_FUNCTION',@getTransformedDistribution_timeSensitive_v1_1,...
                 'USING_METAGEN', 1,...
                 'METAGEN_MAX_TIME', 1000,...
                 'SAVE_NUMBER', key.TIME_LIMIT + 1,...
                 'TIMES_TO_STORE', data.t,...
                 'PERCENT_COMPLETE_INCREMENT', 0.05,...
                 'QUIET_MODE', 0,...
                 'CONV',conv,...
                 'CARRYING_CAPACITY',cc,...
                 'USING_DEATH_RAMP',use_ramp,...
                 'DEATH_RAMP_CENTER',ramp_center,...
                 'DEATH_RAMP_WIDTH',ramp_width,...
                 'DEATH_RAMP_CREATOR',@getSigmoidalRampFunction_v1_0);

% get fitter
fitter = getFitter_v1_4();
fitter_specs = struct('chainer',fitter.options.chainer_PowLawGrowthMeta,...
     'weights',fitter.options.generate_weights_noWeights_lowFilter(1));
% fitter_specs = struct('plotter',fitter.options.plotter_getLinkedPlotter,...
%     'chainer',fitter.options.chainer_PowLawGrowthMeta,...
%     'weights',fitter.options.weights_noWeights);
fitter = getFitter_v1_4(fitter_specs);

% run model
mg_res = conmod_v2_3(key, changes);
mg_res = prePlotCalcs_v1_2(mg_res);

% time comparison points
[mg_fit,~] = fitter.objective(mg_res,data);
axis([data.a(1) max([conv.v2a(mg_res.x(mg_res.max_index)) data.a(end)]) 0 1.1*ceil(max(data.cdf(end,:)))]);
% [~,mg_fit] = fitter.gradient(mg_res,data);
disp(mg_fit);
% plotCDF_v1_3(mg_res,data,[data.a(1) max([conv.v2a(mg_res.x(mg_res.max_index)) data.a(end)]) 0 1.1*ceil(max(data.cdf(end,:)))]);
% mg_fit = -1;