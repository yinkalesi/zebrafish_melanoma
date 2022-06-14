%% getSelectedOptima_Meta_v2_0
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 6/2/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  2.0: using distance_weights and meta-birth fitting
%  2.1: allow submission of key arguments (to help implement changes in
%  ORIGIN_TIME_CALCULATOR)

function [opt_pars,details] = getSelectedOptima_v2_1(data_file_name,name,...
    start_pars,selected,par_low,par_high,conv,times_used,use_ramp,run_time,...
    low_thresh,stage_number,usePara,args)

if(~exist('stage_number','var'))
    stage_number = 3;
end

if(~exist('usePara','var'))
    usePara = true;
end

if(~exist('args','var'))
    args = struct();
end

start_time = clock();
run_time_hrs = run_time;
max_iterations = round(1000*run_time_hrs);

par_names = {'Growth_Parameter1','Growth_Parameter2','Death_Parameter1',...
    'Death_Parameter2','Shed_Parameter1','Meta_Parameter1',...
    'Meta_Parameter2','Carrying_Capacity','Ramp_Center','Ramp_Width'};

% initial guess
kg = start_pars(1:2);
kr = start_pars(3:4);
ks = start_pars(5);
km = start_pars(6:7);
cc = start_pars(8);
ramp_center = start_pars(9);
ramp_width = start_pars(10);

% trunc function used to 'flatten' parameters
tr = @(x) log10(x+1); % log function
tri = @(y) 10.^y-1; % inverse function
% function to translate gradient of normal parameters to flattened
% parameters
trid = @(y) log(10)*10.^y; % help with derivatives
tridp = @(y) log(10)./(1-10.^y).*ones(size(y));

% tr = @(x) log10(x); % log function
% tri = @(y) 10.^y; % inverse function
% % function to translate gradient of normal parameters to flattened
% % parameters
% trid = @(y) log(10)*10.^y; % help with derivatives
% tridp = @(y) log(10)*ones(size(y));

x0 = tr(start_pars(selected));
low =  tr(par_low(selected));
high = tr(par_high(selected));

disp(['init: ' num2str(tri(x0))]);
disp(['low : ' num2str(tri(low))]);
disp(['high: ' num2str(tri(high))]);

data = getExperimentalData_v3_0(data_file_name,conv,times_used);
if(stage_number==1)
    nx = 500;
    dt = 0.05;
elseif(stage_number==2)
    nx = 500;
    dt = 0.05;
else
    nx = 500;
    dt = 0.02;
end
opt_key = getDataKey_v5_1(data,nx,dt,name,start_pars,'uniform');
data = opt_key.inferred_data;

changes = struct('SOLUTION_FUNCTION',@getTransformedDistribution_CTC_v1_0,...
    'USING_METAGEN', 0,...
    'SAVE_NUMBER', opt_key.TIME_LIMIT + 1,...
    'TIMES_TO_STORE', data.t,...
    'PERCENT_COMPLETE_INCREMENT', 0.25,...
    'QUIET_MODE', 1,...
    'CONV',conv,...
    'CARRYING_CAPACITY',cc,...
    'USING_DEATH_RAMP',use_ramp,...
    'DEATH_RAMP_CENTER',ramp_center,...
    'DEATH_RAMP_WIDTH',ramp_width,...
    'DEATH_RAMP_CREATOR',@getSigmoidalRampFunction_v1_0);
args_f = fieldnames(args);
for f = 1:length(args_f)
    field = args_f{f};
    changes.(field) = args.(field);
end

opt_key = changeKey_v3_1(opt_key, changes);

% get fitter
fitter = getFitter_v4_0();
if(stage_number==1)
    if(usePara)
        fitter_specs = struct('plotter',fitter.options.plotter_getNullPlotter,...
            'weights',fitter.options.generate_weights_distanceWeights_TimeDependentLowFilter(low_thresh));
    else
        fitter_specs = struct('plotter',fitter.options.plotter_getLinkedPlotter,...
            'weights',fitter.options.generate_weights_distanceWeights_TimeDependentLowFilter(low_thresh));
        if(evalin('base','exist(''fitter_figure'',''var'')'))
            evalin('base','fitter_figure.delete');
        end
    end
elseif(stage_number==2)
    if(usePara)
        fitter_specs = struct('plotter',fitter.options.plotter_getNullPlotter,...
            'residuals',fitter.options.residuals_meta3,'chainer',fitter.options.chainer_blank,'weights',fitter.options.weights_noWeights);
    else
        fitter_specs = struct('plotter',fitter.options.plotter_getLinkedPlotter,...
            'residuals',fitter.options.residuals_meta3,'chainer',fitter.options.chainer_blank,'weights',fitter.options.weights_noWeights);
        if(evalin('base','exist(''fitter_figure'',''var'')'))
            evalin('base','fitter_figure.delete');
        end
    end
else
    if(usePara)
        fitter_specs = struct('plotter',fitter.options.plotter_getNullPlotter,...
            'weights',fitter.options.generate_weights_distanceWeights_lowFilter(low_thresh));
    else
        fitter_specs = struct('plotter',fitter.options.plotter_getLinkedPlotter,...
            'weights',fitter.options.generate_weights_distanceWeights_lowFilter(low_thresh));
        if(evalin('base','exist(''fitter_figure'',''var'')'))
            evalin('base','fitter_figure.delete');
        end
    end
end
fitter = getFitter_v4_0(fitter_specs);

% need function to generate full parameter set from supplied parameters
if(start_pars(6) ~= 0)
    shed_meta_ratio = start_pars(5)/start_pars(6);
    smr_low = par_low(5)/par_high(6);
    smr_high = par_high(5)/par_low(6);
else
    shed_meta_ratio = 1e3;
    smr_low = 1e2;
    smr_high = 1e4;
end
p0 = tr([start_pars(1:4) shed_meta_ratio start_pars(6:10)]);
plow = tr([par_low(1:4) smr_low par_low(6:10)]);
phigh = tr([par_high(1:4) smr_high par_high(6:10)]);
xin = data.x(data.selector{1}(end));
xfin = data.x(data.selector{end}(end));
% kx = (xfin*log(xfin)-xin*log(xin))/(xfin-xin)-1;
% kx = 0.5*(log(xin)+log(xfin)); 
if(stage_number==1)
    kx = 13;
elseif(stage_number==2)
    kx = 12;
else
    kx = 13.5;
end

pef = @(an,bn,b0) tri(an)*exp(-kx*(tri(bn)-tri(b0))); % scales pre-exponetial factor
pefd1 = @(an,bn,b0) tridp(an)*pef(an,bn,b0); % partial derivative of pef
pefd2 = @(an,bn,b0) -kx*trid(bn)*pef(an,bn,b0); % partial derivative of pef
in2par = @(p) [pef(p(1),p(2),p0(2)) tri(p(2)) ...
    pef(p(3),p(4),p0(4)) tri(p(4)) tri(p(5))*pef(p(6),p(7),p0(7)) ...
    pef(p(6),p(7),p0(7)) tri(p(7:10))];
in2pard = @(p) [pefd1(p(1),p(2),p0(2)) pefd2(p(1),p(2),p0(2)) zeros(1,length(p)-2);...
    zeros(1,1) trid(p(2)) zeros(1,length(p)-2);...
    zeros(1,2) pefd1(p(3),p(4),p0(4)) pefd2(p(3),p(4),p0(4)) zeros(1,length(p)-4);...
    zeros(1,3) trid(p(4)) zeros(1,length(p)-4);...
    zeros(1,4) trid(p(5))*pef(p(6),p(7),p0(7)) tri(p(5))*pefd1(p(6),p(7),p0(7)) tri(p(5))*pefd2(p(6),p(7),p0(7)) zeros(1,length(p)-7);...
    zeros(1,5) pefd1(p(6),p(7),p0(7)) pefd2(p(6),p(7),p0(7)) zeros(1,length(p)-7);...
    zeros(1,6) trid(p(7)) zeros(1,length(p)-7);...
    zeros(1,7) trid(p(8)) zeros(1,length(p)-8);...
    zeros(1,8) trid(p(9)) zeros(1,length(p)-9);...
    zeros(1,9) trid(p(10))];

if(run_time > 0)
    % anonymous function for optimization
    opt_func = @(x) optimizer_v10_0_nograd(in2par(select2in(x,p0,selected)),...
        opt_key,fitter,in2pard(select2in(x,p0,selected)));
    
    brem_opts = struct('del',0.10,'maxIter',max_iterations,'TimeLimit',...
        run_time_hrs*3600,'Tolerance',1e-8,'ToleranceLength',100);
    if(usePara)
        [opt_res,opt_val,flag,opt_out] = stagedOptimization_v1_3(@bremermann_parallel_v2_0,...
            opt_func,p0(selected),plow(selected),phigh(selected),[2 3],brem_opts);
    else
        [opt_res,opt_val,flag,opt_out] = stagedOptimization_v1_3(@bremermann_v2_0,...
            opt_func,p0(selected),plow(selected),phigh(selected),[2 3],brem_opts);
    end
else
    opt_res = x0;
    flag = -2;
    opt_out = struct();
    mg_res = conmod_v3_1(opt_key, changes);
    mg_res = prePlotCalcs_v1_3(mg_res);
    % time comparison points
    [opt_val,~] = fitter.objective(mg_res,data);
end

disp(opt_res);
tres = in2par(select2in(opt_res,p0,selected));
% tres = tri(opt_res);
% tres(1) = tres(1)*exp(-tconst*(tres(2)-gam));
% disp(tres);
fprintf('%f\n',tres(selected));
end_time = clock();

duration_mins = etime(end_time,start_time)/60;
disp([mfilename ' took ' num2str(duration_mins) ' minutes']);

% need simulation log to log results, procedures
res_log = [clock duration_mins opt_key.NUMBER_OF_SIZE_INTERVALS opt_key.TIME_STEP start_pars tres];
dlmwrite(['optimizer_log_' name '.txt'],res_log,'-append');

% convert inputs from optimizer to actual parameters
inputs = opt_out.inputs;
inputs2 = zeros(size(inputs,1),length(tres));
for i = 1:size(inputs,1)
    inputs2(i,:) = in2par(select2in(inputs(i,:),p0,selected));
end
opt_out.inputs = inputs2;

details = struct('opt_out',opt_out,'opt_flag',flag,'opt_val',opt_val,...
    'selected',selected,'optimized_parameters',{par_names(selected)});
opt_pars = tres;
end

function p = select2in(s,p0,selected)
% need to add optimized parameters back into full parameter matrix and
% return the result
p = p0;
p(selected) = s;
end
