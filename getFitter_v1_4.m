function [fitter] = getFitter_v1_4(new_specs)
%% getFitter_v1_3
%  Version 1.3
%  Author: Adeyinka Lesi
%  Date: 3/19/18
%  Project: Tumor Growth, Logarithmic Continuum Form

% return a structure containing fitting functions
%% Version History
%  1.1: using multiplyByColumn for backwards compatibility with Matlab 2015
%  1.2: changing distanceWeights (avoid inf) and adding new weight func
%  1.3: add a null_plotter for when you don't want to plot anything
%  (usefull for parallel runs)
%  1.4: replacing INITIAL_TUMOR_DIST with TIME_ZERO_DIST
%  7/7/18: changed how data cdf is calculated and plotted - it affected the
%  dif used for fitting so previous fits may be imperfect; started using 
%  getRefinedDif, which calculates the dif differentently in the case of
%  primary tumors vs metastasis

fitter = struct();
fit_nosqrt = @(dif,w) w*(dif.^2)/length(dif);
fit_nosqrt_deriv = @(dif,w) 2*w.*dif'/length(dif);
fit_sqrt = @(dif,w) sqrt(w*(dif.^2)/length(dif));
fit_sqrt_deriv = @(dif,w) w.*dif'/fit_sqrt(dif,w)/length(dif);
fit_sqrt2 = @(dif,w) sqrt(w*(dif.^2));
fit_sqrt2_deriv = @(dif,w) w.*dif'/fit_sqrt2(dif,w);

fitter.defaults = struct(...
    'residuals',@getRefinedDif,...
    'weights',@noWeights,...
    'function',fit_sqrt,...
    'derivative',fit_sqrt_deriv,...
    'chainer',@getPowLawGrowthChain,...
    'plot',@semilogx,...
    'plotter',@getPlotter,...
    'function_helper',@getFit);

fitter.options = struct(...
    'residuals_midPointDif',@getDif,...
    'weights_noWeights',@noWeights,...
    'weights_distanceWeights',@distanceWeights,...
    'weights_relativeDistanceWeights',@relativeDistanceWeights,...
    'weights_relativeDistanceWeights_lowPenalty',@relativeDistanceWeights_lowPenalty,...
    'weights_sqrtSizeDistanceWeights_lowPenalty',@sqrtSizeDistanceWeights_lowPenalty,...
    'generate_weights_noWeights_lowFilter',@noWeights_lowFilter_generator,...
    'plotter_getLinkedPlotter',@getLinkedPlotter,...
    'plotter_getNullPlotter',@getNullPlotter,...
    'chainer_PowLawGrowthMeta',@getPowLawGrowthMetaChain,...
    'function_helper_single_tumors',@getFit_withMoments,...
    'fit_function_no_sqrt',fit_nosqrt,...
    'fit_function_no_sqrt_derivative',fit_nosqrt_deriv,...
    'fit_function_sqrt2',fit_sqrt2,...
    'fit_function_sqrt2_derivative',fit_sqrt2_deriv);

specs = fitter.defaults;

if(exist('new_specs','var'))
    changes = fieldnames(new_specs);
    for f = 1:length(changes)
        if(isfield(specs,changes{f}))
            specs.(changes{f}) = new_specs.(changes{f});
        else
            warning(['Unknown specification ' changes{f}]);
        end
    end
end

fitter.specs = specs;
fitter.objective = @(res,data) specs.function_helper(res,data,specs);
fitter.gradient = @(res,data) getGradient(res,data,specs);
end

function [fit,grad] = getFit(res,data,s)
t = res.t2;
% disp(['Times: ' num2str(t)]);
if(any(data.t~=t))
    error('Time variable mismatch during fitting');
end
fit = zeros(1,length(t));
plotter = s.plotter();
plotter.start(res,data,s);
wdifs = cell(1,length(t));
derivs = cell(1,length(t));
xcomps = cell(1,length(t));
for i = 1:length(t)
    [dif,xcomp] = s.residuals(res,data,i);
    w = s.weights(res,data,i,xcomp);
    if(~isempty(dif))
        fit(i) = s.function(dif,w);
        derivs{i} = s.derivative(dif,w);
        wdifs{i} = w.*dif'/length(w);
        xcomps{i} = xcomp;
    else
        % need to use moments to calculate a fit
        % calculate 0th and 1st moment of model
        data_x_low = data.x(data.selector{i}(1));
        i_x_low = res.key.FIELD.indexOf(data_x_low*0.75,res.key);
        rxp = res.xp(i_x_low:end);
        rdx = res.x(i_x_low+1:end)-res.x(i_x_low:end-1);
        dist = res.dist2(i_x_low:end,i);
        rm0 = rdx*dist;
        rm1 = (rdx.*(rxp-0.5))*dist;
        % calculate moments of data
        em0 = data.cdf(end,i);
        em1 = data.x'*data.dist(:,i);
        % calculate an error
        if(rm0 == 0)
            E = 0;
        else
            E = 2*(rm1*em0/rm0/em1-1);
        end
        fit(i) = s.function(E,ones(size(E)));
        derivs{i} = s.derivative(E,ones(size(E)));
        wdifs{i} = E;
        xcomps{i} = em1/em0;
    end
end
plotter.continue(s.plot,res,data,wdifs,xcomps);
plotter.end(res,data,fit);
grad = getFitGradient(res,data,s,derivs,xcomps);
if(any(isnan(fit)))
    warning('Fit gives NaN');
    disp(res.key.PARAMETERS);
end
end

function [fit,grad] = getFit_withMoments(res,data,s)
t = res.t2;
if(any(data.t~=t))
    error('Time variable mismatch during fitting');
end
fit = zeros(1,length(t));
plotter = s.plotter();
plotter.start(res,data,s);
difs = cell(1,length(t));
derivs = cell(1,length(t));
xcomps = cell(1,length(t));
for i = 1:length(t)
    [dif,xcomp] = s.residuals(res,data,i);
    if(~isempty(dif))
        w = s.weights(res,data,i,xcomp);
        fit(i) = s.function(dif,w);
        difs{i} = dif;
        derivs{i} = s.derivative(dif,w);
        xcomps{i} = xcomp;
    else
        % calculate 0th and 1st moment of model
        rxp = res.xp;
        rdx = res.x(2:end)-res.x(1:end-1);
        dist = res.dist2(:,i);
        rm0 = rdx*dist;
        rm1 = (rdx.*(rxp-0.5))*dist;
        % calculate moments of data
        em0 = data.cdf(end,i);
        em1 = data.x'*data.dist(:,i);
        % calculate an error
        if(rm0 == 0)
            E = 0; % during metagen time...
        else
            E = rm1/rm0-em1/em0;
        end
        fit(i) = E.^2;
        derivs{i} = 2*E;
        difs{i} = E;
        xcomps{i} = data.x;
    end
end
plotter.continue(s.plot,res,data,difs,xcomps);
plotter.end(res,data,fit);
grad = getFitGradient(res,data,s,derivs,xcomps);
end

function [grad,fit] = getGradient(res,data,s)
t = res.t2;
if(any(data.t~=t))
    error('Time variable mismatch during fitting');
end
fit = zeros(1,length(t));
plotter = s.plotter();
plotter.start(res,data,s);
difs = cell(1,length(t));
derivs = cell(1,length(t));
xcomps = cell(1,length(t));
for i = 1:length(t)
    [dif,xcomp] = s.residuals(res,data,i);
    w = s.weights(res,data,i,xcomp);
    fit(i) = s.function(dif,w);
    difs{i} = dif;
    derivs{i} = s.derivative(dif,w);
    xcomps{i} = xcomp;
end
[grad,chains] = getFitGradient(res,data,s,derivs,xcomps);
plotter.chain(s.plot,res,data,chains,xcomps);
plotter.chainEnd(res,data,grad,chains);
end

function [grad,chains] = getFitGradient(res,data,s,derivs,xcomps)
[chains] = s.chainer(res,data,xcomps);
grad = zeros(length(chains),length(data.t));
for i = 1:length(data.t)
    for j = 1:length(chains)
        grad(j,i) = -derivs{i}*chains{j}{i};
    end
end
end

function [dif,xcomp] = getDif(res,data,i)
% data ccdf
dx = data.x(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
if(isfield(data,'cdf_adj'))
    p_adj = data.cdf_adj(i);
else
    p_adj = 0;
end
% dccdf = dcdf(end)+p_adj - dcdf;
dccdf1 = p_adj+dcdf(end)-[0; dcdf(1:end-1)];

% model ccdf
imin = 1;
max_dxi = res.key.FIELD.indexOf(dx(end),res.key);
imax = max(res.max_index,max_dxi);
mcdf = res.dist_cum2(imin:imax+1,i);
mccdf = mcdf(end) - mcdf;

xcomp = 0.5*(dx(2:end)+dx(1:end-1));
inds = res.key.FIELD.indexOf(xcomp,res.key) - imin + 1;

dif = 0.5*(dccdf1(1:end-1)+dccdf1(2:end)) - mccdf(inds);
end

function [dif,xcomp] = getRefinedDif(res,data,i)
% data ccdf
dx = data.x(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
if(isfield(data,'cdf_adj'))
    p_adj = data.cdf_adj(i);
else
    p_adj = 0;
end
% dccdf = dcdf(end)+p_adj - dcdf;
dccdf1 = p_adj+dcdf(end)-[0; dcdf(1:end-1)];

% model ccdf
imin = 1;
max_dxi = res.key.FIELD.indexOf(dx(end),res.key);
imax = max(res.max_index,max_dxi);
mcdf = res.dist_cum2(imin:imax+1,i);
mccdf = mcdf(end) - mcdf;

xcomp = 0.5*(dx(2:end)+dx(1:end-1));
inds = res.key.FIELD.indexOf(xcomp,res.key) - imin + 1;

plev = data.cdf(data.selector{1}(end),1);
iplev = find(dccdf1<=plev,1,'first');
if(isempty(iplev))
    iplev = length(dccdf1);
end
dif = zeros(length(inds),1);
dif(1:iplev-1) = 0.5*(dccdf1(1:iplev-1)+dccdf1(2:iplev)) - mccdf(inds(1:iplev-1));
dif(iplev:end) = dccdf1(iplev+1:end) - mccdf(inds(iplev:end));
end

function [w] = noWeights(res,data,i,xcomp)
w = ones(1,length(xcomp))/length(xcomp);
end

function [w] = distanceWeights(res,data,i,xcomp)
da = data.a(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
del = da(2:end)-da(1:end-1);
dis_ave = (dcdf(2:end)-dcdf(1:end-1))./del;
w = 1./dis_ave;
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end

function [w] = relativeDistanceWeights(res,data,i,xcomp)
da = data.a(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
del = da(2:end)-da(1:end-1);
dis_ave = (dcdf(2:end)-dcdf(1:end-1))./del;
denom = dis_ave.*xcomp;
w = 1./denom.^2;
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end

function [w] = relativeDistanceWeights_lowPenalty(res,data,i,xcomp)
da = data.a(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
del = da(2:end)-da(1:end-1);
dis_ave = (dcdf(2:end)-dcdf(1:end-1))./del;
denom = dis_ave.*(xcomp+1000*data.lowa./xcomp); % penalizes really small and really large tumors
w = 1./denom.^2;
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end

function [w] = sqrtSizeDistanceWeights_lowPenalty(res,data,i,xcomp)
da = data.a(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
del = da(2:end)-da(1:end-1);
dis_ave = (dcdf(2:end)-dcdf(1:end-1))./del;
denom = dis_ave.*sqrt(xcomp+1000*data.lowa./xcomp); % penalizes really small and really large tumors
w = 1./denom.^2;
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end

function [w] = noWeights_lowFilter(res,data,i,xcomp,limit)
da = data.a(data.selector{i});
da_mid = 0.5*(da(2:end)+da(1:end-1));
w = ones(size(da_mid));
w(da(1:end-1)<limit) = 0;

if(~isempty(w))
    w(end) = 1; % make sure at least one nonzero val
end
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end

function [weight_func] = noWeights_lowFilter_generator(limit)
    weight_func = @(res,data,i,xcomp) noWeights_lowFilter(res,data,i,xcomp,res.key.CONV.v2a(limit));
end

function [plotter] = getPlotter()
    plotter = struct();
    plotter.start = @(res,data,s) figure;
    plotter.continue = @plotterContinue2;
    plotter.chain = @plotterChain;
    plotter.end = @plotterEnd;
    plotter.chainEnd = @plotterChainEnd;
end

function [] = plotterContinue(plot_func,res,data,i,dif,xcomp)
if(any(i==data.plot_indices))
    da = data.a(data.selector{i});
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(i);
    else
        p_adj = 0;
    end
%     dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
    dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
    imin = 1;
    imax = res.max_index;
    ma = res.key.CONV.v2a(res.x(imin:imax+1));
    mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
    colors = hsv(length(data.t));
    dac = zeros(2*length(da),1);
    dac(1:2:end-1) = da;
    dac(2:2:end) = da;
    dccdfc = zeros(2*length(da),1);
    dccdfc(1:2:end-1) = dccdf1;
    dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
    plot_func(dac,dccdfc,'o:',ma,mccdf,'-',res.key.CONV.v2a(xcomp),abs(dif),'.',...
        'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
        'LineWidth',1.5,'MarkerSize',4);
    hold on;
end
end

function [] = plotterContinue2(plot_func,res,data,difs,xcomps)
for i = 1:length(data.t)
    dif = difs{i};
    xcomp = xcomps{i};
    if(any(i==data.plot_indices))
        da = data.a(data.selector{i});
        if(isfield(data,'cdf_adj'))
            p_adj = data.cdf_adj(i);
        else
            p_adj = 0;
        end
%         dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
        dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
        imin = 1;
        imax = res.max_index;
        ma = res.key.CONV.v2a(res.x(imin:imax+1));
        mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
        colors = hsv(length(data.t));
            dac = zeros(2*length(da),1);
    dac(1:2:end-1) = da;
    dac(2:2:end) = da;
    dccdfc = zeros(2*length(da),1);
    dccdfc(1:2:end-1) = dccdf1;
    dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
    plot_func(dac,dccdfc,'o:',ma,mccdf,'-',res.key.CONV.v2a(xcomp),abs(dif),'.',...
        'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
        'LineWidth',1.5,'MarkerSize',4);
    hold on;
    end
end
end

function [] = plotterChain(plot_func,res,data,chains,xcomps)
for i = 1:length(data.t)
    xcomp = xcomps{i};
    if(any(i==data.plot_indices))
        da = data.a(data.selector{i});
        if(isfield(data,'cdf_adj'))
            p_adj = data.cdf_adj(i);
        else
            p_adj = 0;
        end
%         dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
        dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
        imin = 1;
        imax = res.max_index;
        ma = res.key.CONV.v2a(res.x(imin:imax+1));
        mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
        colors = hsv(length(data.t));
        dac = zeros(2*length(da),1);
        dac(1:2:end-1) = da;
        dac(2:2:end) = da;
        dccdfc = zeros(2*length(da),1);
        dccdfc(1:2:end-1) = dccdf1;
        dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
        plot_func(dac,dccdfc,'o:',ma,mccdf,'-',...
            'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
            'LineWidth',1.5,'MarkerSize',4);
        hold on;
        for j = 1:length(chains)
            chain = abs(chains{j}{i});
            ylim = max([data.cdf(end,:) res.dist_cum2(end,:)]);
            ylim = 10*floor(0.1*ylim)+10;  
            chain = chain/max(chain)*ylim;
            plot_func(res.key.CONV.v2a(xcomp),chain,'--',...
                'Color',max(0,(0.6-0.05*j))*colors(i,:));
        end
    end
end
end

function [] = plotterEnd(res,data,fit)
hold off;
title(res.key.TITLE);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
% get legend
leg = cell(1,3*length(data.plot_indices));
for i = 1:length(data.plot_indices)
    leg{3*i-2} = ['Data, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{3*i-1} = ['Model, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{3*i} = ['fit = ' sprintf('%1.3f',fit(data.plot_indices(i)))];
end
legend(leg,'Location','Best');
alow = data.lowa;
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.a > alow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
ylim = max([data.cdf(end,:)-data.cdf(ia1,:) res.dist_cum2(end,:)-res.dist_cum2(ix1,:)]);
ylim = floor(ylim)+1; 
axis([alow inf 0 ylim]);
end

function [] = plotterChainEnd(res,data,grad,chains)
hold off;
title(res.key.TITLE);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
% get legend
nl = length(chains)+2;
leg = cell(1,nl*length(data.plot_indices));
for i = 1:length(data.plot_indices)
    leg{nl*i-(nl-1)} = ['Data, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{nl*i-(nl-2)} = ['Model, t = ' num2str(data.t(data.plot_indices(i)))];
    for j = 3:nl
        leg{nl*i-(nl-j)} = ['chain' num2str(j-2) ' = ' num2str(grad(j-2,i))];
    end
end
legend(leg,'Location','Best');
alow = 1;
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.a > alow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
ylim = max([data.cdf(end,:)-data.cdf(ia1,:) res.dist_cum2(end,:)-res.dist_cum2(ix1,:)]);
ylim = floor(ylim)+1; 
axis([alow inf 0 ylim]);
end

function [plotter] = getLinkedPlotter()
    plotter = struct();
    plotter.start = @linkedPlotterStart;
    plotter.continue = @linkedPlotterContinue;
    plotter.end = @linkedPlotterEnd;
    plotter.chain = @linkedPlotterContinue;
    plotter.chainEnd = @(res,data,grad,chains) linkedPlotterEnd(res,data,grad);
end

function [] = linkedPlotterStart(res,data,s)
if(evalin('base','exist(''fitter_figure'',''var'')'))
    if(~evalin('base','fitter_figure.isvalid'))
        evalin('base','fitter_figure = figure;');
        initializeLinkedPlot(s.plot,res,data);
    end
else
    evalin('base','fitter_figure = figure;');
    initializeLinkedPlot(s.plot,res,data);
end
end

function [] = initializeLinkedPlot(plot_func,res,data)
% initialize plot
evalin('base','fitter_figure.Position = [1 1 950 680];');
ra = res.key.CONV.v2a(res.x(1:end));
assignin('base','fitter_ra',ra);
assignin('base','fitter_res',addByColumn(...
    -res.dist_cum2(1:end,:),res.dist_cum2(end,:)));
dat_a = cell(1,length(data.t));
dat_cdf = cell(1,length(data.t));
for i = 1:length(dat_cdf)
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(i);
    else
        p_adj = 0;
    end
    dat_ai = data.a(data.selector{i});
    dat_ai2 = zeros(2*length(dat_ai),1);
    dat_ai2(1:2:end-1) = dat_ai;
    dat_ai2(2:2:end) = dat_ai;
    dat_a{i} = dat_ai2;
    dcdfi = p_adj+data.cdf(data.selector{i}(end),i)-...
        [0; data.cdf(data.selector{i}(1:end-1),i)];
    dcdfi2 = zeros(2*length(dat_ai),1);
    dcdfi2(1:2:end-1) = dcdfi;
    dcdfi2(2:2:end) = [dcdfi(2:end); p_adj];
    dat_cdf{i} = dcdfi2;
end
assignin('base','fitter_da',dat_a);
assignin('base','fitter_dat',dat_cdf);

% make plots
colors = hsv(length(data.t));
plotcol = [colors;0.6*colors];
assignin('base','fitter_colors',plotcol);
evalin('base','figure(fitter_figure);');
evalin('base','plot(0,0);');
evalin('base','set(fitter_figure.CurrentAxes,''ColorOrder'',fitter_colors);');
evalin('base','set(fitter_figure.CurrentAxes,''NextPlot'',''replacechildren'');');
evalin('base',[func2str(plot_func) ...
    '(fitter_ra,fitter_res,''LineWidth'',1.25);']);
hold on;
for i = 1:length(data.t)
    evalin('base',[func2str(plot_func) '(fitter_da{' num2str(i) '},'...
        'fitter_dat{' num2str(i) '},''o--'',''LineWidth'',0.75);']);
end
hold off;
title(res.key.TITLE);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
evalin('base','linkdata(fitter_figure)');

% store fit variable
evalin('base','fitter_min_fit = inf;');
end

function [] = linkedPlotterContinue(plot_func,res,data,difs,xcomps)
ra = res.key.CONV.v2a(res.x(1:end));
assignin('base','fitter_ra',ra);
assignin('base','fitter_res',addByColumn(...
    -res.dist_cum2(1:end,:),res.dist_cum2(end,:)));
end

function [] = linkedPlotterEnd(res,data,fit)
evalin('base','figure(fitter_figure);');
title(res.key.TITLE);
% get legend
leg = cell(1,2*length(data.t));
for i = 1:length(data.t)
    leg{i+length(data.t)} = ['Data, t = ' num2str(data.t(i))];
    leg{i} = ['Model Fit =' sprintf(' %1.3f',fit(:,i))];
end
legend(leg,'Location','Best');
alow = max(data.lowa,data.a(1));
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.a > alow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
ylim_res = max(res.dist_cum2(end,:)-res.dist_cum2(ix1,:));
ylim_dat = max(data.cdf(end,:)-data.cdf(ia1,:));
ylim = max(ylim_dat,min(ylim_dat*1.5,ylim_res));
ylim = floor(ylim)+1; 
axis([alow inf 0 ylim]);
% axis([1e0 2e4 0 6]);

evalin('base','refreshdata(fitter_figure);');
pause(0.1);

% update saved fit
if(evalin('base',['fitter_min_fit > ' num2str(max(fit))]))
    assignin('base','fitter_min_fit',max(fit));
    assignin('base','fitter_min_parameters',res.key.PARAMETERS);
end
evalin('base','disp([''MINIMUM FIT: '' num2str(fitter_min_fit)]);');

end

function [plotter] = getNullPlotter()
    plotter = struct();
    plotter.start = @(res,data,s) fprintf('Null Plot ');
    plotter.continue = @(plot_func,res,data,difs,xcomps) fprintf('-');
    plotter.end = @(res,data,fit) fprintf(' Ended\n');
    plotter.chain = @(plot_func,res,data,difs,xcomps) fprintf('-');
    plotter.chainEnd = @(res,data,grad,chains) fprintf('-\n');
end

function [chains] = getPowLawGrowthChain(res,data,xcomps)
t = data.t;
dt = 0.5*res.key.TIME_STEP;
chain1 = cell(1,length(t));
chain2 = cell(1,length(t));
x_init = res.key.TIME_ZERO_SIZES;
traj = getSizeTrajectory_v2_1(res.key,x_init,dt,t(end));
ti = round(t/dt)+1;
xt = traj(:,ti);
for i = 1:length(t)
    xcomp = xcomps{i};
    chain1{i} = zeros(size(xcomp));
    chain2{i} = zeros(size(xcomp));
    ixn = res.key.FIELD.indexOf(xcomp,res.key);
    pn = max(0,res.dist2(ixn,i)); % no negatives please (happens for unstable solutions)
    for n = 1:length(xcomp)
        xn = xcomp(n);
        % find tumors closest to xn size (fluctuations in the size of
        % these tumors affects the distribution the most)
        inear = find(xt(:,i)>xn,1);
        % find fluctuation in size due to parameter change
        if(isempty(inear))
            % xn is larger than all xt and only largest tumor matters
            dxda = getParameterDeriv1(xt(end,i),x_init(end),t(i),res.key);
            dxdb = getParameterDeriv2(xt(end,i),x_init(end),t(i),res.key);
        elseif(inear == 1)
            % xn is smaller that all xt and only smallest tumor and
            % metastases matter
%             ntraj = getSizeTrajectory_v2_1(res.key,xn,-dt,t(i));
            c1 = res.key.GROWTH_PARAMETER1;
            c2 = res.key.GROWTH_PARAMETER2;
            t_max = (xn^(1-c2)-1)/c1/(1-c2);
            if(t(i)>t_max)
                dxda = getParameterDeriv1(xn,1,t_max,res.key);
                dxdb = getParameterDeriv2(xn,1,t_max,res.key);
            else
                xn0 = (xn^(1-c2)-c1*(1-c2)*t(i)).^(1/(1-c2));
                dxda = getParameterDeriv1(xn,xn0,t(i),res.key);
                dxdb = getParameterDeriv2(xn,xn0,t(i),res.key);
            end
        else
            % xn in sandwiched in between two closest tumors in xt
            dxdas = getParameterDeriv1(xt(inear-1:inear,i),...
                x_init(inear-1:inear),t(i),res.key);
            dxdbs = getParameterDeriv2(xt(inear-1:inear,i),...
                x_init(inear-1:inear),t(i),res.key);
            wei = abs(1./(xt(inear-1:inear,i)-xn));
            wei = wei/sum(wei);
            dxda = wei'*dxdas;
            dxdb = wei'*dxdbs;
        end
        % find fluctuation in ccdf
        chain1{i}(n) = pn(n)*dxda;
        chain2{i}(n) = pn(n)*dxdb;
    end
end
chains = {chain1,chain2};
end

function [chains] = getPowLawGrowthMetaChain(res,data,xcomps)
% parameters: [growth1, growth2, death1, death2, meta1, meta2]
pnums = [1 2 5 6];
t = data.t;
dt = 0.5*res.key.TIME_STEP;
chain1 = cell(1,length(t));
chain2 = cell(1,length(t));
chain3 = cell(1,length(t));
chain4 = cell(1,length(t));
x_init = res.key.TIME_ZERO_SIZES;
p_init = res.key.TIME_ZERO_DIST;
traj = getSizeTrajectory_v2_1(res.key,x_init,dt,t(end)); % trajectories for primaries
[traj1,t1] = getSizeTrajectory_v2_1(res.key,1,dt,t(end)); % trajectory from size 1
% summation parameters for traj1
dx1 = traj1(2:end)-traj1(1:end-1);
x1 = traj1(1:end-1)+0.5*dx1;
vx1 = getV(x1,res.key);
vtraj1 = getV(traj1,res.key);
% parameter derivatives for velocity
dv1dp = zeros(length(pnums),length(x1));
dvt1dp = zeros(length(pnums),length(traj1));
dDt1dp = zeros(length(pnums),length(traj1));
for pn = 1:length(pnums)
    dv1dp(pn,:) = getV_PowLawParamDeriv(x1,pnums(pn),res.key);
    dvt1dp(pn,:) = getV_PowLawParamDeriv(traj1,pnums(pn),res.key);
    dDt1dp(pn,:) = getD_PowLawParamDeriv(traj1,pnums(pn),res.key);
end
% calculate del_t function
ddtdp = multiplyByColumn(dv1dp,-(dx1./vx1.^2));
dtdp = zeros(length(pnums),length(x1)+1);
for di = 1:length(ddtdp)
    dtdp(:,di+1) = dtdp(:,di)+ddtdp(:,di);
end
% calculate dx1dp/v function (change in trajectory from size 1 divided by velocity)
ddx1dp_v = multiplyByColumn(dv1dp,dt./vx1);
dx1dp_v = zeros(length(pnums),length(x1)+1);
for di = 1:length(ddx1dp_v)
    dx1dp_v(:,di+1) = dx1dp_v(:,di)+ddx1dp_v(:,di);
end
% calculate dxjdp/v function (change in trajectory of primaries divided by
% velocity)
np = size(traj,1); % number of primary tumors
nx = size(traj,2); % number of time points/trajectory points
dxj = traj(:,2:end)-traj(:,1:end-1);
xj = traj(:,1:end-1)+0.5*dxj;
vxj = zeros(np,nx-1);
vtraj = zeros(np,nx);
% calculate parameter derivative summation terms
dvdpj = zeros(length(pnums),nx-1,np);
ddxjdp_v = zeros(length(pnums),nx-1,np);
for ni = 1:np
    % parameter derivatives for velocity
    for pn = 1:length(pnums)
        dvdpj(pn,:,ni) = getV_PowLawParamDeriv(xj(ni,:),pnums(pn),res.key);
    end
    vxj(ni,:) = getV(xj(ni,:),res.key);
    vtraj(ni,:) = getV(traj(ni,:),res.key);
    ddxjdp_v(:,:,ni) = multiplyByColumn(dvdpj(:,:,ni),dt./vxj(ni,:));
end
% tabulate parameter derivatives for each primary
dxjdp_v = zeros(length(pnums),nx,np);
for di = 1:nx-1
    for ni = 1:np
        dxjdp_v(:,di+1,ni) = dxjdp_v(:,di,ni)+ddxjdp_v(:,di,ni);
    end
end

% factor for change in distribution due to change in metastasis generation
% from primaries
kmp = zeros(np,nx);
for ni = 1:np
    kmp(ni,:) = res.key.RATES.meta(traj(ni,:));
end
sum_km = sum(kmp);
dkmpdp = zeros(length(pnums),nx,np);
for pn = 1:length(pnums)
    for ni = 1:np
        pam1 = getLKM_PowLawParamDeriv(traj(ni,:),pnums(pn),res.key);
        pam2 = res.key.META_PARAMETER2*dxjdp_v(pn,:,ni).*vtraj(ni,:)./traj(ni,:);
        dkmpdp(pn,:,ni) = kmp(ni,:).*(pam1+pam2);
    end
end
dLMdp = zeros(length(pnums),nx);
for pn = 1:length(pnums)
    dkmpdp_sum = zeros(1,nx);
    for ni = 1:np
        dkmpdp_sum = dkmpdp_sum + dkmpdp(pn,:,ni);
    end
    dlnkg1dp = getLKG1_PowLawParamDeriv(traj(ni,:),pnums(pn),res.key);
    dLMdp(pn,:) = dkmpdp_sum./sum_km - dlnkg1dp;
end

% factor for change in how velocity changes affect movement of metastasis
dvdx = getV_deriv(traj1,res.key); % derivative of velocity in terms of x
dLVdp = zeros(length(pnums),nx);
for pn = 1:length(pnums)
    % need variable dvt1dp instead of dv1dp because dvt1dp is calculated at
    % traj1 while dv1dp is calculated at x1 (intermediate points)
    dLVdp(pn,:) = dvt1dp(pn,:)./vtraj1+dvdx.*dx1dp_v(pn,:);
end

% factor for change in attenuation due to diffusivity (modeling this
% attenuation as 0.5/sqrt(D(x)*t), which is only a very rough approximation
% since D changes with position - but the argument that this is good enough
% is that we only care about the relative change in attenuation so even if
% the value of the estimate is wrong, the derivative of the log may be more
% reasonable
Dtraj1 = getD(traj1,res.key);
dDdx = getD_deriv(traj1,res.key);
dLAdp = zeros(length(pnums),nx);
for pn = 1:length(pnums)
    dLAdp(pn,:) = -0.5*(dDt1dp(pn,:)+dDdx.*vtraj1.*dx1dp_v(pn,:))./Dtraj1;
end

ti = round(t/dt)+1;
xt = traj(:,ti);
itraj1 = res.key.FIELD.indexOf(traj1,res.key);
itraj1 = min(itraj1,size(res.dist2,1)); % to avoid exceeding index size
for i = 1:length(t)
    xcomp = xcomps{i};
    chain1{i} = zeros(size(xcomp));
    chain2{i} = zeros(size(xcomp));
    chain3{i} = zeros(size(xcomp));
    chain4{i} = zeros(size(xcomp));
    ixn = res.key.FIELD.indexOf(xcomp,res.key);
    pn = max(0,res.dist2(ixn,i)); % no negatives please (happens for unstable solutions)
    % velocity at xcomp
    vn = getV(xcomp,res.key);
    Dn = getD(xcomp,res.key);
    % values for quadrature
    % t1 -> times (edge values)
    % x1 -> trajectories (midpoint values)
    % vx1 -> velocities (midpoint values)
    % ix1 -> indices for x1 in distribution
    ptraj1 = max(0,res.dist2(itraj1,i));
    
    for n = 1:length(xcomp)
        xn = xcomp(n);
        % 1) find change due to movement of primaries
        % find tumors closest to xn size (fluctuations in the size of
        % these tumors affects the distribution the most)
        inear = find(xt(:,i)>xn,1);
        % find fluctuation in size due to parameter change
        if(isempty(inear))
            % xn is larger than all xt and only largest tumor matters
            dxda = vtraj(end,ti(i))*dxjdp_v(1,ti(i),end);
            dxdb = vtraj(end,ti(i))*dxjdp_v(2,ti(i),end);
            % no metastasis effect
            dtdpn = zeros(length(pnums),1);
            dMdpn = zeros(length(pnums),1);
        elseif(inear ~= 1)
            % xn in sandwiched in between two closest tumors in xt
            x1 = xt(inear-1,i);
            x2 = xt(inear,i);
            dxdas = vtraj(inear-1:inear,ti(i)).*[dxjdp_v(1,ti(i),inear-1) dxjdp_v(1,ti(i),inear)]';
            dxdbs = vtraj(inear-1:inear,ti(i)).*[dxjdp_v(2,ti(i),inear-1) dxjdp_v(2,ti(i),inear)]';
            %             test1 = getParameterDeriv1(xt(inear-1:inear,i),...
            %                 x_init(inear-1:inear),t(i),res.key);
            %             test2 = getParameterDeriv2(xt(inear-1:inear,i),...
            %                 x_init(inear-1:inear),t(i),res.key);
            %             disp([dxdas test1 dxdbs test2]);
            % need to account for number of tumors at each position as well
            % as distance from xn to interpolate the proper derivative
            % estimates
            p1 = p_init(inear-1);
            p2 = p_init(inear);
            if(p1==0 || p2==0)
                error('Unexpected zero value');
            end
            f1 = (x2-xn)/p2;
            f2 = (xn-x1)/p1;
            wei = [f1 f2]/(f1+f2);
            dxda = wei*dxdas;
            dxdb = wei*dxdbs;
            % no metastasis effect
            dtdpn = zeros(length(pnums),1);
            dMdpn = zeros(length(pnums),1);
        else
            % xn is smaller that all xt and we need to incorporate
            % metastasis generation effects
            
            % given posibilty that xn is just smaller than the size of the
            % smallest tumor, still want to take affect of smallest tumors
            % into account but at the same time, this effect needs to fade
            % out for situations where xn is a lot smaller and smallest
            % tumor has no effect (will assume a Gaussian probability
            % desity for the smallest tumor)
            attenuation = exp(-(xn-xt(1,i))^2/4/Dn(n)); % not sure if this value is the right one (since no time dependence reflected here)
            dxda = vn(n)*dxjdp_v(1,ti(i),1)*attenuation;
            dxdb = vn(n)*dxjdp_v(2,ti(i),1)*attenuation;
            
            % 2) find change due to movement of metastasis wave
            % calculate growth time to get to xcomp(n)
            x1n_i = find(xn<traj1,1); % t1(x1n_i) is larger than the time it takes to grow to xn from 1
            if(isempty(x1n_i))
                % a metastasis hasn't had the time to grow to this size
                dtdpn = zeros(length(pnums),1);
                dMdpn = zeros(length(pnums),1);
            elseif(x1n_i==1)
                % this should not happen since traj1(1) == 1 and xn>=1
                disp('Converter predicts data point indicates a tumor of a size lower than 1');
                dtdpn = dtdp(:,1); % initial values are all zeros
                dMdpn = zeros(length(pnums),1);
            else
                % will use linear interpolation
                dtdpn = dtdp(:,x1n_i-1)+(dtdp(:,x1n_i)-dtdp(:,x1n_i-1))*...
                    (xn-traj1(x1n_i-1))/(traj1(x1n_i)-traj1(x1n_i-1));
                
                % 3) find change due to differences in metastasis generation
                % time to grow form size 1 to xn (linear interpolation)
                dtn = t1(x1n_i-1)+(t1(x1n_i)-t1(x1n_i-1))*...
                    (xn-traj1(x1n_i-1))/(traj1(x1n_i)-traj1(x1n_i-1));
                % need to perform  quadrature: tsys=system refrence time,
                % tgrow=growth time=t-tsys are important variables
                idtn = round(dtn/dt)+1; % index for low index of integration for tgrow
                it_dtn = ti(i)-idtn+1; % index for high limit of integration for tsys
                if(it_dtn>0)
                    % tgrow in integrated in reverse => fliplr(dtn:t)
                    itgrow = fliplr(idtn:ti(i));
                    % velocity of tgrow=fliplr(dtn:t)
                    vgrow = vtraj1(itgrow);
                    % probability densities for tgrow
                    pgrow = ptraj1(itgrow)';
                    
                    % a) relative change in metastasis generation of primaries
                    % for each parameter -> dMdp
                    % b) relative change in velocity -> dLVdpL
                    % c) relative change in diffusive attenuation -> dLAdp
                    
                    % add up the changes
                    dMdpn = dt*(vgrow.*pgrow)*...
                        (dLMdp(:,1:it_dtn)+dLVdp(:,itgrow)+dLAdp(:,itgrow))';
                else
                    % not enough time has passed for metastasis to grow
                    % this large or very small growth time - can aproximate 
                    % the derivative as zero
                    dMdpn = zeros(length(pnums),1);
                end
            end
        end
        
        % find fluctuation in ccdf
        chain1{i}(n) = pn(n)*dxda-pn(n)*dtdpn(1)*vn(n)+dMdpn(1);
        chain2{i}(n) = pn(n)*dxdb-pn(n)*dtdpn(2)*vn(n)+dMdpn(2);
        chain3{i}(n) = -pn(n)*dtdpn(3)*vn(n)+dMdpn(3);
        chain4{i}(n) = -pn(n)*dtdpn(4)*vn(n)+dMdpn(4);
    end
end
chains = {chain1,chain2,chain3,chain4};
end

function [dxda] = getParameterDeriv1(xt,x0,t,key)
b = key.GROWTH_PARAMETER2;
dxda = t*xt.^b;
end

function [dxdb] = getParameterDeriv2(xt,x0,t,key)
a = key.GROWTH_PARAMETER1;
b = key.GROWTH_PARAMETER2;
term1 = xt.*log(xt);
term2 = -a*t*xt.^b;
term3 = -(xt.^b).*(x0.^(1-b)).*log(x0);
dxdb = (term1+term2+term3)/(1-b);
end

function [v] = getV(x,key)
    v = key.RATES.growth(x)-key.RATES.death(x)-key.RATES.shed(x);
end

function [dvdx] = getV_deriv(x,key)
    alp = key.GROWTH_PARAMETER2;
    bet = key.DEATH_PARAMETER2;
    gam = key.META_PARAMETER2;
    dvdx = 0.5*(alp*key.RATES.growth(x)-bet*key.RATES.death(x)-gam*key.RATES.shed(x))./x;
end

function [dvda] = getV_PowLawParamDeriv(x,pnum,key)
    switch pnum
        case 1
            a = key.GROWTH_PARAMETER1;
            if(a>0)
                dvda = key.RATES.growth(x)/a;
            else
                dvda = zeros(size(x));
            end
        case 2
            dvda = key.RATES.growth(x).*log(x);
        case 3
            b = key.DEATH_PARAMETER1;
            if(b>0)
                dvda = -key.RATES.death(x)/b;
            else
                dvda = zeros(size(x));
            end
        case 4
            dvda = -key.RATES.death(x).*log(x);
        case 5
            c = key.META_PARAMETER1;
            if(c>0)
                dvda = -key.RATES.shed(x)/c;
            else
                dvda = zeros(size(x));
            end
        case 6
            dvda = -key.RATES.shed(x).*log(x);
    end
end

function [D] = getD(x,key)
    D = 0.5*(key.RATES.growth(x)+key.RATES.death(x)+key.RATES.shed(x));
end

function [dDdx] = getD_deriv(x,key)
    alp = key.GROWTH_PARAMETER2;
    bet = key.DEATH_PARAMETER2;
    gam = key.META_PARAMETER2;
    dDdx = 0.5*(alp*key.RATES.growth(x)+bet*key.RATES.death(x)+gam*key.RATES.shed(x))./x;
end

function [dDda] = getD_PowLawParamDeriv(x,pnum,key)
    switch pnum
        case 1
            a = key.GROWTH_PARAMETER1;
            if(a>0)
                dDda = 0.5*key.RATES.growth(x)/a;
            else
                dDda = zeros(size(x));
            end
        case 2
            dDda = 0.5*key.RATES.growth(x).*log(x);
        case 3
            b = key.DEATH_PARAMETER1;
            if(b>0)
                dDda = 0.5*key.RATES.death(x)/b;
            else
                dDda = zeros(size(x));
            end
        case 4
            dDda = 0.5*key.RATES.death(x).*log(x);
        case 5
            c = key.META_PARAMETER1;
            if(c>0)
                dDda = 0.5*key.RATES.shed(x)/c;
            else
                dDda = zeros(size(x));
            end
        case 6
            dDda = 0.5*key.RATES.shed(x).*log(x);
    end
end

function [dMda_km] = getLKM_PowLawParamDeriv(x,pnum,key)
    % note: this is calculating derivative divided by km(x), i.e. the
    % derivative of the log
    switch pnum
        case 1
            dMda_km = zeros(size(x));
        case 2
            dMda_km = zeros(size(x));
        case 3
            dMda_km = zeros(size(x));
        case 4
            dMda_km = zeros(size(x));
        case 5
            c = key.META_PARAMETER1;
            if(c>0)
                dMda_km = ones(size(x))/c;
            else
                dMda_km = zeros(size(x));
            end
        case 6
            dMda_km = log(x);
    end
end

function [dkg1da_kg1] = getLKG1_PowLawParamDeriv(x,pnum,key)
    % note: this is calculating derivative divided by km(x), i.e. the
    % derivative of the log
    switch pnum
        case 1
            a = key.GROWTH_PARAMETER1;
            dkg1da_kg1 = ones(size(x))/a;
        case 2
            dkg1da_kg1 = zeros(size(x));
        case 3
            dkg1da_kg1 = zeros(size(x));
        case 4
            dkg1da_kg1 = zeros(size(x));
        case 5
            dkg1da_kg1 = zeros(size(x));
        case 6
            dkg1da_kg1 = zeros(size(x));
    end
end