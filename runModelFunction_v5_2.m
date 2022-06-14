%% runModelFunction_v5_2
%  Version 5.2
%  Author: Adeyinka Lesi
%  Date: 3/22/21
%  Project: Tumor Growth, Logarithmic Continuum Form

%% Version History
%  2.0: using getDataKey_v4_0 (new way to calculate time_zero_sizes)
%  3.0: allows working with no data using constructData_v1_0
%  4.0: implementing getTransformedDistribution_CTC_v1_0
%  4.1: changeable grid size and time step
%  4.2: using getFitter_v3_1 (with MMD fit func)
%  4.3: using getFitter_v3_2 (with Adjs MMD)
%  4.4: using LLS
%  4.5: using getFitter_v3_3 (generate_weights_noWeights_TimeDependentLowFilter)
%  4.8: distanceWeights
%  4.9: mirroring getSelectedOptima_v1_16 with stage_number
%  5.0: implementing time zero size depth > 1 and using getDataKey_v5_1
%  12/31/20: changed to noWeights for stage3 to align with
%  getSelectedOptima_v2_3
%  5.1: using weights_noWeights_MetaFilter
%  5.2: adding data extrapolation (3/22/21), switching to firstIsZero mode

function [mg_res,data,mg_fit] = runModelFunction_v5_2(datafile,params,conv,Nx,dt,select,args,stage_number,extp_mat)

% if datafile is a string with '=', we construct a data file using the
% times sizes given => t1,t2=x11,x12;x21,x22
ieq = find('='==datafile,1,'first');
isc = find(';'==datafile,1,'first');
% set name
% name by combinining data file with time
if(isempty(ieq))
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
else
    % no data case
    clockdat = clock;
    if(isempty(isc))
        isc = length(datafile)+1;
    end
    size_str = datafile(ieq+1:isc-1);
    size_str(','==size_str) = '_';
    name = ['x0_eq_' size_str '_' sprintf('%02i%02i%4i-%02i%02i',clockdat([2 3 1 4 5]))];
end

if(~exist('conv','var'))
    conv = getNullConverter();
end
if(~exist('select','var'))
    select = [];
end
if(~exist('args','var'))
    args = struct();
end

% load data
if(isempty(ieq))
    if(~isempty(select))
        data = getExperimentalData_v3_2(datafile,conv,select);
        %         data = getExperimentalData_v3_3(datafile,conv,select,extp_mat);
    else
        data = getExperimentalData_v3_2(datafile,conv);
        %         data = getExperimentalData_v3_3(datafile,conv,[],extp_mat);
    end
else
    times = str2num(datafile(1:ieq-1));
    sizes = str2num(datafile(ieq+1:end));
    dist = zeros(length(sizes(:)),length(times));
    for i = 1:length(times)
        dist((i-1)*size(sizes,2)+1:i*size(sizes,2),i) = 1;
    end
    vec = sizes';
    vec = vec(:);
    data = constructData_v1_0(times,vec,dist,conv);
end

% for later use
cc = params(8);
ramp_center = params(9);
ramp_width = params(10);
use_ramp = ramp_center > 0;

% key = getKey_metagen_v1_1([79.616798681379152   0.569437688675373],[0.182851433668764   0.947026744394068],ks,[ks*1e-3 gamma], ...
%     data.t(end),data.t(end)*10,data.x,data.dist(:,1),2000,primary_init,data.lowx,name);

key = getDataKey_v5_3(data,Nx,dt,name,params,'uniform'); % using getKey_backtrack_v5_0
data = key.inferred_data;

changes = struct('SOLUTION_FUNCTION',@getTransformedDistribution_CTC_v1_1,...
                 'USING_METAGEN', 0,...
                 'SAVE_NUMBER', min(5*key.TIME_LIMIT + 1,101),...
                 'TIMES_TO_STORE', data.t,...
                 'PERCENT_COMPLETE_INCREMENT', 0.05,...
                 'QUIET_MODE', 0,...
                 'CONV',conv,...
                 'CARRYING_CAPACITY',cc,...
                 'USING_DEATH_RAMP',use_ramp,...
                 'DEATH_RAMP_CENTER',ramp_center,...
                 'DEATH_RAMP_WIDTH',ramp_width,...
                 'DEATH_RAMP_CREATOR',@getSigmoidalRampFunction_v1_0);
changed = fieldnames(args);
if(~isempty(changed))
    for f = 1:length(changed)
        field = changed{f};
        changes.(field) = args.(field);
    end
end

% % get fitter
% fitter = getFitter_v3_2();
% % not using weights by setting limit to zero: fitter.options.generate_weights_MMDWeights(conv.a2v(6))
% fitter_specs = struct('weights',fitter.options.generate_weights_MMDWeights(0),...
%         'residuals',fitter.options.generate_residuals_AdjMMD(datafile,select),...
%         'function',fitter.options.fit_function_mmd);
% fitter = getFitter_v3_2(fitter_specs);

% get fitter
low_thresh = 0;
fitter = getFitter_v4_1();
if(stage_number==1)
    fitter_specs = struct('weights',fitter.options.generate_weights_distanceWeights_TimeDependentLowFilter(low_thresh));
elseif(stage_number==2)
    %     fitter_specs = struct('weights',fitter.options.generate_weights_noWeights_lowFilter(low_thresh));
    fitter_specs = struct('residuals',fitter.options.residuals_meta5,...
        'chainer',fitter.options.chainer_blank,'weights',fitter.options.weights_noWeights);
else
    fitter_specs = struct('weights',fitter.options.generate_weights_noWeights_lowFilter(low_thresh));
    %     fitter_specs = struct('weights',fitter.options.generate_weights_distanceWeights_lowFilter(low_thresh));
    %     fitter_specs = struct('plotter',fitter.options.plotter_getDifPlotter,...
    %         'weights',fitter.options.weights_countWeights);
    
    %     data_file_name = datafile;
    %     times_used = select;
    %     fitter_specs = struct('weights',fitter.options.generate_weights_MMDWeights(low_thresh),...
    %         'residuals',fitter.options.generate_residuals_AdjMMD(data_file_name,times_used),...
    %         'function',fitter.options.fit_function_mmd);
end
fitter = getFitter_v4_1(fitter_specs);


% run model
mg_res = conmod_v3_1(key, changes);
mg_res = prePlotCalcs_v1_2(mg_res);
data = mg_res.key.inferred_data;
mg_res.data = data;
mg_res.fitter = fitter;

% % plots
% figure;
% subplot(2,1,1)
% k0pk1 = (mg_res.key.CTC_PARAMETER1+mg_res.key.CTC_PARAMETER2);
% plot(mg_res.t3,mg_res.ctc,mg_res.t3,mg_res.ctc1,mg_res.t3,mg_res.num_shed/k0pk1,'--');
% legend('TM1','RK4','C_{ss}');
% xlabel('Time');
% ylabel('CTC Concentration');
% subplot(2,1,2)
% plot(mg_res.t3,mg_res.num_meta/mg_res.key.FIELD.dt);
% xlabel('Time');
% ylabel('Metastasis Generation');
% % time comparison points
[mg_fit,~] = fitter.objective(mg_res,data);
% [mg_grad,mg_fit] = fitter.gradient(mg_res,data);
disp(mg_fit);
% ylim1 = 0;
% ylim2 = max(mg_res.dist_cum(end),1.1*ceil(max(data.cdf(end,:))));
% loc1_t0 = find(mg_res.dist_cum(:,1)>0.05,1,'first');
% loc1_tf = find(mg_res.dist_cum(:,end)>0.05,1,'first');
% xlim1 = min([data.a(1) mg_res.x(loc1_t0) mg_res.x(loc1_tf)]);
% loc2_t0 = find(mg_res.dist_cum(end,1)-mg_res.dist_cum(:,1)<0.05,1,'first');
% loc2_tf = find(mg_res.dist_cum(end,end)-mg_res.dist_cum(:,end)<0.05,1,'first');
% xlim2 = max([data.a(end) mg_res.x(loc2_t0) mg_res.x(loc2_tf)]);
% axis([xlim1 xlim2 ylim1 ylim2]);