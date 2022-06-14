%% getFitter_v4_1
%  Version 4.1
%  Author: Adeyinka Lesi
%  Date: 1/21/21
%  Project: Tumor Growth, Logarithmic Continuum Form
% return a structure containing fitting functions
%% Version History
%  1.1: using multiplyByColumn for backwards compatibility with Matlab 2015
%  1.2: changing distanceWeights (avoid inf) and adding new weight func
%  1.3: add a null_plotter for when you don't want to plot anything
%  (usefull for parallel runs)
%  1.4: replacing INITIAL_TUMOR_DIST with TIME_ZERO_DIST
%  2.0: implementing robust, accurate gradient calculator, further changes
%  to dif
%  3.0: works with getTransformedDistribution_CTC_v1_0
%  3.1: implementing alternate objective function (MMD like function); also
%  switched default to fit_sqrt2 (the weights already average dif^2 so you
%  don't need to dived by length(dif) to average)
%  12/18/19: there was an error in getRefinedDif where the ccdf didn't
%  necessarily start at infinity - this messed up the fitting at times.
%  3.2: implementing advanced version of MMD that accounts for fish being
%  taken out of the distribution
%  3.3: implementing noWeights_lowThresh with variable cutoff depending on
%  timepoint. 3/19/20: modified plotter to have dif stand out more
%  4.0: using inferred data in fit
%  4.1: adding primary fit option; changing plot lower y limit to 0.01 to
%  accomodate loglog plots, adding new meta fit options, adding plotting
%  options

function [fitter] = getFitter_v4_1(new_specs)

fitter = struct();
fit_nosqrt = @(dif,w) w*(dif.^2)/length(dif);
fit_nosqrt_deriv = @(dif,w) 2*w.*dif'/length(dif);
fit_sqrt = @(dif,w) sqrt(w*(dif.^2)/length(dif));
fit_sqrt2 = @(dif,w) sqrt(w*(dif.^2));
fit_mmd = @(dif,w) w*abs(dif);
% made expression more complicated to account for case where dif = 0 - now
% give proper value instead of NaN for that case
fit_sqrt_mod = @(dif,w) any(dif)*fit_sqrt(dif,w)+1-any(dif);
fit_wdif_mod = @(dif,w) any(dif)*w.*dif'+sqrt(w)*(1-any(dif));
fit_sqrt_deriv = @(dif,w) fit_wdif_mod(dif,w)/fit_sqrt_mod(dif,w)/length(dif);
fit_sqrt2_mod = @(dif,w) any(dif)*fit_sqrt2(dif,w)+1-any(dif);
fit_sqrt2_deriv = @(dif,w) fit_wdif_mod(dif,w)/fit_sqrt2_mod(dif,w);

fitter.defaults = struct(...
    'residuals',@getRefinedDif,...
    'weights',@noWeights,...
    'function',fit_sqrt2,...
    'derivative',fit_sqrt2_deriv,...
    'chainer',@getFullChain,...
    'plot',@semilogx,...
    'plotter',@getPlotter,...
    'function_helper',@getFit);

fitter.options = struct(...
    'residuals_midPointDif',@getDif,...
    'residuals_mmd',@getMMDDif,...
    'residuals_meta',@getMetastasisDif_byOriginTime,...
    'residuals_meta2',@getMetastasisDif_byMetCount,...
    'residuals_meta3',@getMetastasisDif_byCountAndSize,...
    'residuals_meta4',@getMetastasisDif_byCountSizeJump,...
    'residuals_meta5',@getMetastasisDif_byCountAndSize2,...
    'residuals_meta6',@getMetastasisDif_byCountSizeJump2,...
    'residuals_meta7',@getMetastasisDif_byCountSizeJump3,...
    'residuals_meta8',@getMetastasisDif_byCountSizeJump4,...
    'residuals_meta9',@getMetastasisDif_byCountAndSize3,...
    'generate_residuals_AdjMMD',@getAdjMMDDif_generator,...
    'weights_noWeights',@noWeights,...
    'weights_distanceWeights',@distanceWeights,...
    'weights_pdfWeights',@pdfWeights,...
    'weights_countWeights',@countWeights,...
    'weights_noWeights_PrimaryFilter',@noWeights_PrimaryFilter,...
    'weights_noWeights_MetaFilter',@noWeights_MetaFilter,...
    'weights_distanceWeights_PrimaryFilter',@distanceWeights_PrimaryFilter,...
    'weights_pdfWeights_PrimaryFilter',@pdfWeights_PrimaryFilter,...
    'generate_weights_distanceWeights_lowFilter',@distanceWeights_lowFilter_generator,...
    'generate_weights_pdfWeights_lowFilter',@pdfWeights_lowFilter_generator,...
    'generate_weights_distanceWeights_TimeDependentLowFilter',@distanceWeights_TimeDependentLowFilter_generator,...
    'weights_relativeDistanceWeights',@relativeDistanceWeights,...
    'weights_relativeDistanceWeights_lowPenalty',@relativeDistanceWeights_lowPenalty,...
    'weights_sqrtSizeDistanceWeights_lowPenalty',@sqrtSizeDistanceWeights_lowPenalty,...
    'wieghts_MMDWeights',@getMMDWeights,...
    'wieghts_heightAdjusted',@heightAdjustedWeights,...
    'generate_weights_noWeights_lowFilter',@noWeights_lowFilter_generator,...
    'generate_weights_heightAdjustedWeights_lowFilter',@heightAdjustedWeights_lowFilter_generator,...
    'generate_weights_noWeights_TimeDependentLowFilter',@noWeights_TimeDependentLowFilter_generator,...
    'generate_weights_heightAdjustedWeights_TimeDependentLowFilter',@heightAdjustedWeights_TimeDependentLowFilter_generator,...
    'generate_weights_MMDWeights',@getMMDWeights_generator,...
    'generate_weights_PDFMMDWeights',@getPDFMMDWeights_generator,...
    'plotter_getLinkedPlotter',@getLinkedPlotter,...
    'plotter_getNullPlotter',@getNullPlotter,...
    'plotter_getDifPlotter',@getDifPlotter,...
    'plotter_noInference',@getPlotterNoInference,...
    'function_helper_single_tumors',@getFit_withMoments,...
    'fit_function_no_sqrt',fit_nosqrt,...
    'fit_function_no_sqrt_derivative',fit_nosqrt_deriv,...
    'fit_function_sqrt_per_point',fit_sqrt,...
    'fit_function_sqrt__per_point_derivative',fit_sqrt_deriv,...
    'fit_function_sqrt2',fit_sqrt2,...
    'fit_function_sqrt2_derivative',fit_sqrt2_deriv,...
    'fit_function_mmd',fit_mmd,...
    'chainer_blank',@getBlankChain);

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
        wdifs{i} = w.*dif.^2';
        xcomps{i} = xcomp;
    elseif(~isempty(data.selector{i}))
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
    else
        fit(i) = 0;
        derivs{i} = 0;
        wdifs{i} = 0;
        xcomps{i} = 1;
    end
end
plotter.continue(s.plot,res,data,wdifs,xcomps);
plotter.end(res,data,fit);
grad = getFitGradient(res,data,s,derivs,xcomps);
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

function [grad,fit,chains] = getGradient(res,data,s)
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
        grad(j,i) = -derivs{i}*chains{j,i};
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
dccdf = dcdf(end)+p_adj - dcdf;

% model ccdf
imin = 1;
max_dxi = res.key.FIELD.indexOf(dx(end),res.key);
imax = max(res.max_index,max_dxi);
mcdf = res.dist_cum2(imin:imax+1,i);
mccdf = mcdf(end) - mcdf;

xcomp = 0.5*(dx(2:end)+dx(1:end-1));
inds = res.key.FIELD.indexOf(xcomp,res.key) - imin + 1;
xedge1 = res.x(inds);
xedge2 = res.x(inds+1);
mccdfx = mccdf(inds)+(xcomp-xedge1)./(xedge2-xedge1).*(mccdf(inds+1)-mccdf(inds));

dif = 0.5*(dccdf(1:end-1)+dccdf(2:end)) - mccdfx;
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

if(isempty(dcdf))
    dif = [];
    xcomp = [];
else
    % dccdf = dcdf(end)+p_adj - dcdf;
    dccdf1 = p_adj+dcdf(end)-[0; dcdf(1:end-1)];
    
    % model ccdf
    imin = 1;
    % max_dxi = res.key.FIELD.indexOf(dx(end),res.key);
    imax = res.max_index; % if this value goes out of bounds, look at HARD_TUMOR_SIZE_LIMIT compared to sizes in data
    mcdf = res.dist_cum2(imin:imax+1,i);
    mccdf = mcdf(end) - mcdf;
    
    xcomp = 0.5*(dx(2:end)+dx(1:end-1));
    inds = min(res.max_index,res.key.FIELD.indexOf(xcomp,res.key)) - imin + 1;
    xedge1 = res.x(inds)';
    xedge2 = res.x(inds+1)';
    mccdfx = mccdf(inds)+(xcomp-xedge1)./(xedge2-xedge1).*(mccdf(inds+1)-mccdf(inds));
    
    plev = data.cdf(data.selector{1}(end),1);
    iplev = find(dccdf1<=plev,1,'first');
    if(isempty(iplev))
        iplev = length(dccdf1);
    end
    dif = zeros(length(inds),1);
    dif(1:iplev-1) = 0.5*(dccdf1(1:iplev-1)+dccdf1(2:iplev)) - mccdfx(1:iplev-1);
    dif(iplev:end) = dccdf1(iplev+1:end) - mccdfx(iplev:end);
end
end

function obj = getAdjMMDDif_generator(datfiles,select)
    obj = @(res,data,i) getAdjMMDDif(res,data,i,datfiles,select);
end

function [dif,xcomp] = getAdjMMDDif(res,data,i,datfiles,select)
% the idea: represent the discrete distribution with Gaussian functions and
% use the difference in number density (integrated over space) to calculate
% the objective function
% note 1: integrating the abs(difference) -> Levy distance
% note 2: this method has issues when the number of fish used to calculate
% the distribution changes over time. This adjusted version checks the
% trajectory file to see which fish are still present

% model dist
dist1 = abs(res.dist2);
t2 = res.t2;
it3 = round(t2/(res.t3(2)-res.t3(1))+1);
x = res.x';
xcomp = res.xp';
dxs = res.key.FIELD.dxs';
m_x0 = res.key.TIME_ZERO_SIZES;
m_dist0 = res.key.TIME_ZERO_DIST;
cramp = res.complementary_ramp(it3);
% data dist
datx = max(1,data.x(data.selector{i})');
datdist = data.dist(data.selector{i},i);
datdist(datx<=1) = 0;
datcdf = data.cdf(data.selector{i},i);

% construct approximate data distribution using Gaussian kernel
cg = @(x1,xi,Di,w) (1/2-erf((x1-xi)./sqrt(2*Di))/2)*w;
gm = @(x1,xi,Di) -diff(cg(x1,xi,Di,datdist))./diff(x1);
gm2 = @(x1,xi,Di) -diff(cg(x1,xi,Di,ones(length(xi),1)))./diff(x1);

idatx = min(length(dxs),res.key.FIELD.indexOf(datx,res.key));
dx_datx = dxs(idatx)';
Di_corr = dx_datx.*2;
Di = Di_corr.*t2(i).*max(0.5*(res.key.RATES.growth(datx)+(1-cramp(i))*res.key.RATES.death(datx)+res.key.RATES.shed(datx)),1);

gmx1 = gm(x,datx,Di)+gm(x,2-datx,Di);

% construct adjusted model distribution
countfile = [datfiles(1:strfind(datfiles,'area')-1) 'counts.txt'];
counts = dlmread(countfile);
sizefile = [datfiles(1:strfind(datfiles,'area')-1) 'sizes.txt'];
sizes = dlmread(sizefile);
% disp(['   Fish Counts: ' num2str(fish_counts)]);
% disp(['Primary Counts: ' num2str(primary_counts)]);
% disp(['  Tumor Counts: ' num2str(tumor_counts)]);
% disp(['  Cured Counts: ' num2str(total_cured_counts)]);
% disp([' Active Counts: ' num2str(active_counts)]);
% disp(['Active Primary: ' num2str(nprims_active)]);
norm_mark = datfiles(strfind(datfiles,'over_')+5:find(datfiles=='.',1,'last')-1);
norm_selector = struct('count',1,'fish',2,'prim',3,'tum',4,'active',5,'actprim',6);
renorm1 = counts(norm_selector.(norm_mark),select(i));
renorm2 = counts(norm_selector.(norm_mark),select(1));
% it is best to use data normalized by active primary with MMD - since the
% fitting makes tumor by tumor comparisons, it makes the most sense when
% the height (ccdf value) we get to the last primary at each time point 
% stays constant, which happens when normalizing by primary. Essentially,
% to make tumor by tumor comparison, I split the model distribution in two
% based on when the density due to primary tumors end so that I can make
% tumor to tumor comparisons with the primary tumor. If the data
% distribution is not normalized by primary, the height at which the split
% occurs changes with time and the ccdf curves don't match (really what
% happens is that you have to adjust the model curve up or down so you get
% a good shape fit but the zero is, of course, now not at zero)
% nprim_factor = counts(6,select(1))/counts(6,select(i))*renorm1/renorm2;

traj = getSizeTrajectory_v3_4(res.key,m_x0,res.key.FIELD.dt,t2(i),res.capacity_factor);
mxprim = traj(:,end)';
% trajind = min(length(dxs),res.key.FIELD.indexOf(mxprim,res.key));
m_count0 = round(sum(m_dist0)*renorm2);
% traj_met = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,t2(i),res.capacity_factor);
% mx_met = traj_met(end);
% metind = min(length(dxs),res.key.FIELD.indexOf(mx_met,res.key));
metind = find(res.dist_cum2(end,i)-res.dist_cum2(:,i)>=m_count0/renorm2,1,'last');
if(isempty(metind))
    metind = 1;
end

% look at sizes file to find out when primaries end up disappearing
init_sizes = sizes(:,select(1));
[~,ord_init] = sort(init_sizes);
prim_start = find(init_sizes(ord_init)>1,1,'first');
prim_sizes = sizes(ord_init(prim_start:end),:);
prim_present = prim_sizes > 0;

cur_sizes_ind = find(prim_present(:,select(i)));

mxprim2 = zeros(1,m_count0);
m_dist2 = round(m_dist0*renorm2);
nmx2 = 0;
for mxi = 1:length(m_dist2)
    mxprim2(nmx2+(1:m_dist2(mxi))) = mxprim(mxi);
    nmx2 = nmx2+m_dist2(mxi);
end

mx_remove_ind = setdiff(1:m_count0,cur_sizes_ind);
mx_remove = mxprim2(mx_remove_ind);
mx_remove_loc = min(res.key.FIELD.indexOf(mx_remove,res.key),length(dxs));

% construct approximate data distribution using Gaussian kernel
dx_mxr = dxs(mx_remove_loc)';
Dmx_corr = dx_mxr.*2;
D_mx = Dmx_corr.*t2(i).*max(0.5*(res.key.RATES.growth(mx_remove)+(1-cramp(i))*res.key.RATES.death(mx_remove)+res.key.RATES.shed(mx_remove)),1);
gm_remove = gm2(x,mx_remove,D_mx);

% dx_tj1 = dxs(trajind(1));
% D1_corr = dx_tj1.*2;
% D1 = D1_corr.*t2(i).*max(0.5*(res.key.RATES.growth(mxprim(1))+(1-cramp(i))*res.key.RATES.death(mxprim(1))+res.key.RATES.shed(mxprim(1))),1);
% spread1 = sqrt(2*D1);
% splitind2 = res.key.FIELD.indexOf(max(mxprim(1)-spread1,1),res.key);
% splitind1 = metind; %round(0.5*(metind+trajind(1)));
% splitind = splitind1; %max(splitind1,splitind2);
splitind = metind;
dist_mod = max(0,[renorm1*dist1(1:splitind-1,i); renorm2*dist1(splitind:end,i)]-gm_remove);
ntum_gt1 = round(sum(datdist(datx>1))*renorm1);
lastmeta = find(round((datcdf(end)-[0;datcdf(1:end-1)])*renorm1)>=ntum_gt1,1,'last');
lastind = idatx(lastmeta);

% truncate model distribution based on number of tumors in data
cdf1 = [0; cumsum(renorm1*gmx1.*dxs)];
ccdf1 = cdf1(end)-cdf1;
cdf2 = [0; cumsum(dist_mod.*dxs)];
ccdf2 = cdf2(end)-cdf2;
ntum_gtli1 = ccdf1(lastind);
if(ccdf2(1)>ntum_gtli1)
    mod_start = find(ccdf2>ntum_gtli1,1,'last');
    dat_start = lastind;
    ntum_comp = ntum_gtli1;
else
    mod_start = 1;
    % adjust dist_mod to reflect missing tumors
    dist_mod(1) = dist_mod(1)+(ntum_gtli1-ccdf2(1))/dxs(1);
    dat_start = lastind;
    ntum_comp = ntum_gtli1;
end
dif = zeros(size(gmx1));
dif(dat_start:end) = renorm1*gmx1(dat_start:end);
dif(mod_start:end) = dif(mod_start:end) - dist_mod(mod_start:end);
dif = dif/ntum_comp;

% construct witness RKHS function using Gaussian kernel
% % dif = (renorm1*gmx1-dist_mod)/ntum_gt1; % normalize based on number of tumors
% % dif(1:lastind-1) = 0; % only compare where tumors are present
% % calculate normalizations
% cdf1 = [0; cumsum(renorm1*gmx1.*dxs)];
% cdf2 = [0; cumsum(dist_mod.*dxs)];
% ntum_gtli1 = cdf1(end)-cdf1(lastind);
% ntum_gtli2 = cdf2(end)-cdf2(lastind);
% ntum_adj = abs(log(ntum_gtli2/ntum_gtli1))/log(2);%abs((ntum_gtli1-ntum_gtli2)/ntum_gtli1);
% uni_dist = 1./(x(end)-x(lastind))*ones(size(dxs));
% if(ntum_gtli2>0)
%     dif = abs(renorm1*gmx1/ntum_gtli1-dist_mod/ntum_gtli2)+ntum_adj*uni_dist; % normalize so integrals from last ind are 1
% else
%     % assume a uniform distribution
%     dif = abs(renorm1*gmx1/ntum_gtli1-uni_dist)+ntum_adj*uni_dist;
% end
% dif(1:lastind-1) = 0; % only compare where tumors are present

% % Plots for troubleshooting
% % recalc cdf
% cdf2 = [0; cumsum(dist_mod.*dxs)];
% ccdf2 = cdf2(end)-cdf2;
% figure;
% loglog(xcomp,renorm1*gmx1/ntum_comp,xcomp,renorm2*dist1(:,i)/ntum_comp,xcomp,dist_mod/ntum_comp,'g',xcomp,abs(dif),'m--');
% title('Rescaled Distribution Comparison');
% legend('Bin-Avg','Model','Modified Model','Dif','Location','Best');
% axis([1 inf 1e-15 inf]);
% hold on; 
% plot(x([dat_start mod_start]),abs(dif([dat_start mod_start])),'*');
% % plot(x([splitind1 splitind2]),abs(dif([splitind1 splitind2])),'o',...
% %     x([dat_start mod_start]),abs(dif([dat_start mod_start])),'*');
% hold off;
% 
% % figure;
% % loglog(xcomp,renorm1/ntum_gtli1*gmx1,xcomp,dist_mod/ntum_gtli2,'g',xcomp,abs(dif),'m--');
% % title('Rescaled Distribution Comparison 2');
% % legend('Bin-Avg','Modified Model','Dif','Location','Best');
% % axis([1 inf 1e-15 inf]);
% 
% % cdf1 = [0; cumsum(gmx1.*dxs)]*renorm1/ntum_gt1;
% % cdf2 = [0; cumsum(dist_mod.*dxs)]/ntum_gt1;
% % ccdf1 = (cdf1(end)-cdf1)/ntum_gtli1;
% % ccdf2 = (cdf2(end)-cdf2)/ntum_gtli2;
% cmmd2 = [0; cumsum(abs(dif).*dxs)];
% ccdf3 = (res.dist_cum2(end,i)-res.dist_cum2(:,i))*renorm1;
% 
% figure;
% semilogx(x,ccdf1/ntum_comp,datx,ccdf1(idatx)/ntum_comp,'o',x,ccdf3/ntum_comp,...
%     x,ccdf2/ntum_comp,x,cmmd2,'--',x(dat_start),ccdf1(dat_start)/ntum_comp,'*',...
%     x(mod_start),ccdf2(mod_start)/ntum_comp,'*',x(splitind),ccdf2(splitind)/ntum_comp,'s');
% title(['CCDF with Obj, t=' num2str(i) ', score=' num2str(round(cmmd2(end),2))]);
% legend('Data','Data Pts','Model','Modified Model','MMD','Location','Best');
% fprintf('Error Score: %3.2f\n',cmmd2(end));
% 
% close; close;%close;
end

function [dif,xcomp] = getMMDDif(res,data,i)
% the idea: represent the discrete distribution with Gaussian functions and
% use the difference in number density (integrated over space) to calculate
% the objective function
% note 1: integrating the abs(difference) -> Levy distance
% note 2: this method has issues when the number of fish used to calculate
% the distribution changes over time. Since we are comparing distributions
% directly, this punishes the model for accurately simulating the growth of
% initial tumors and pushes the parameters to be smaller (esentially a zero
% distribution is better than a wrong distribution)

t2 = res.t2;
x = res.x';
xcomp = res.xp';
dxs = res.key.FIELD.dxs';
datx = max(1,data.x(data.selector{i})');
datdist = data.dist(data.selector{i},i);

% construct approximate data distribution using Gaussian kernel
% g = @(x,xi,Di) (exp(-(x-xi).^2./Di/2)./sqrt(2*pi()*Di))*datdist;
cg = @(x,xi,Di) (1/2-erf((x-xi)./sqrt(2*Di))/2)*datdist;
gm = @(x,xi,Di) -diff(cg(x,xi,Di))./diff(x);

idatx = res.key.FIELD.indexOf(datx,res.key);
dx_datx = dxs(idatx)';
Di_corr = dx_datx.*2;
Di = Di_corr.*t2(i).*max(0.5*(res.key.RATES.growth(datx)+res.key.RATES.death(datx)+res.key.RATES.shed(datx)),1);

gmx1 = gm(x,datx,Di);

% construct witness RKHS function using Gaussian kernel
dif = gmx1-res.dist2(:,i);

% figure;
% loglog(xcomp,gmx1,xcomp,res.dist2(:,i),xcomp,abs(dif),'--');
% title('Distribution Comparison');
% legend('Bin-Avg','Model','Difference','Location','Best');
% axis([1 inf 1e-15 inf]);
end

function [dif,xcomp] = getMetastasisDif_byOriginTime(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);
% calculate dif and xcomp
xcomp = orig_times;
xcompi = round(xcomp/(res.t3(2)-res.t3(1)))+1;
dif = met_cumdist-cum_meta(xcompi)';

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',xcomp,abs(dif),'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byMetCount(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births_i = zeros(1,length(orig_times));
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births_i(j) = locj;
    else
        % set birth time to last time in simulation
        mod_births_i(j) = length(cum_meta);
    end
end

% calculate dif and xcomp
xcomp = met_cumdist;
dif = orig_times-res.t3(mod_births_i)';

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(cumsum(res.num_meta),res.t3,stepcdf,stepx,'o-',xcomp,abs(dif),'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountAndSize(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end);
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';

% calculate xcomp and difs based on birth times and size at detection
xcomp = met_cumdist;
dif1 = (orig_times-mod_births)./mod_births;
dif2 = (det_sizes-est_size_atDet)./det_sizes;
dif = sqrt(dif1.^2+dif2.^2);

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountSizeJump(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results
% 1/29/21:
% this is getMetastasisDif_byCountAndSize with an additional term added in
% to avoid big jumps in origin time

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        % extrapolate birth times for remaining metastases
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end); %set to arbitrary large value to error function evaluate this as a bad parameter set
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';

% calculate xcomp and difs based on birth times and size at detection
xcomp = met_cumdist;
dif1 = (orig_times-mod_births)./mod_births; %since I'm dividing by mod_births, note this emphasizes fitting the earliest metastases
dif2 = (det_sizes-est_size_atDet)./det_sizes;
% dif3 is meant to penalize appararent (large) discontinuities in birth
% time curve. We do this by comparing the slope of the origin time of real
% metastasis to that expected from the model and penalizing large slopes
dif3 = 0.5*[(orig_times(2:end)-orig_times(1:end-1))-(mod_births(2:end)-mod_births(1:end-1)); 0]./[mod_births(2:end)-mod_births(1:end-1); 1]; % multiplying by 0.25 so dif3 doesn't overwhelm dif1
dif = sqrt(dif1.^2+dif2.^2+dif3.^2);
% a note about normalization: I'm dividing my of of the terms in the
% subtraction when calculating the difs so that the different difs end up
% or potentially end up around the same magnitude; they denote a percent
% deviation so it is reasonable to add them together

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountAndSize2(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results
% 2/1/21: changing normalization of dif1 and dif2 since dividing by
% met_count makes it so that fitting low counts are empasized too much

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end);
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
% traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
% est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';
traj_full = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,res.t(end),res.capacity_factor);
ind_det_sizes = sum(traj_full<=det_sizes,2); % this finds the index of first time size becomes larger than det_size for each det_size; this is inefficient but don't want to write a loop
time_to_det_sizes = res.t3(ind_det_sizes)';

% calculate xcomp and difs based on birth times and size at detection
% to combine dif1 and dif2, we need to make sure the differences are in the
% same units/scales - otherwise one might be several magnitudes larger than
% the other. One solution is to convert the size difference represented by
% dif2 to a time difference using growth parameters
xcomp = met_cumdist;
dif1 = (orig_times-mod_births);
dif2 = (time_to_det_sizes-growth_times);
dif = sqrt(dif1.^2+dif2.^2);

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountSizeJump2(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results
% 2/1/21: changing normalization of dif1 and dif2 since dividing by
% met_count makes it so that fitting low counts are empasized too much
% adding dif3 to disinsentivize large jumps

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end);
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
% traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
% est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';
traj_full = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,res.t(end),res.capacity_factor);
ind_det_sizes = sum(traj_full<=det_sizes,2); % this finds the index of first time size becomes larger than det_size for each det_size; this is inefficient but don't want to write a loop
time_to_det_sizes = res.t3(ind_det_sizes)';

% calculate xcomp and difs based on birth times and size at detection
% to combine dif1 and dif2, we need to make sure the differences are in the
% same units/scales - otherwise one might be several magnitudes larger than
% the other. One solution is to convert the size difference represented by
% dif2 to a time difference using growth parameters
xcomp = met_cumdist;
dif1 = (orig_times-mod_births);
dif2 = (time_to_det_sizes-growth_times);
dif3 = [(orig_times(2:end)-orig_times(1:end-1))-(mod_births(2:end)-mod_births(1:end-1)); 0]; 

dif = sqrt(dif1.^2+dif2.^2+dif3.^2);

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountSizeJump3(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results
% 2/1/21: changing normalization of dif1 and dif2 since dividing by
% met_count makes it so that fitting low counts are empasized too much
% adding dif3 to disinsentivize large jumps
% 2/8/21: the last observed time for each metastasis is chosen as the point
% to backtrack from :<

% calculate birth time step function
meta_origin = zeros(size(res.key.META_ORIGIN(:,end)));
for j = 1:length(meta_origin)
    [i1,i2] = find(data.alternate(end).origin_time_loc==j,1,'first');
    if(data.alternate(end).origin_time_corr(i1,end)>=0)
        meta_origin(j) = data.alternate(end).origin_time_corr(i1,end);
    else
        meta_origin(j) = data.alternate(end).origin_time_corr(i1,i2);
    end
end
[~,mord]=sort(meta_origin);
orig_times = meta_origin(mord);
det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
met_cumdist = cumsum(res.key.META_ORIGIN(mord,end-1));
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end);
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
% traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
% est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';
traj_full = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,res.t(end),res.capacity_factor);
ind_det_sizes = sum(traj_full<=det_sizes,2); % this finds the index of first time size becomes larger than det_size for each det_size; this is inefficient but don't want to write a loop
time_to_det_sizes = res.t3(ind_det_sizes)';

% calculate xcomp and difs based on birth times and size at detection
% to combine dif1 and dif2, we need to make sure the differences are in the
% same units/scales - otherwise one might be several magnitudes larger than
% the other. One solution is to convert the size difference represented by
% dif2 to a time difference using growth parameters
xcomp = met_cumdist;
dif1 = (orig_times-mod_births);
dif2 = (time_to_det_sizes-growth_times);
dif3 = [(orig_times(2:end)-orig_times(1:end-1))-(mod_births(2:end)-mod_births(1:end-1)); 0]; 

dif = sqrt(dif1.^2+dif2.^2+dif3.^2);

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountSizeJump4(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results
% 2/1/21: changing normalization of dif1 and dif2 since dividing by
% met_count makes it so that fitting low counts are empasized too much
% adding dif3 to disinsentivize large jumps
% 2/8/21: the last observed time for each metastasis is chosen as the point
% to backtrack from :<
% 2/15/21: using alternate tumor count combining versions 2 and 3

% calculate birth time step function
% count version 1
meta_origin1 = res.key.META_ORIGIN(:,end);
[~,mord1]=sort(meta_origin1);
orig_times1 = meta_origin1(mord1);
met_cumdist1 = cumsum(res.key.META_ORIGIN(mord1,end-1));

% count version 2
meta_origin2 = zeros(size(res.key.META_ORIGIN(:,end)));
for j = 1:length(meta_origin2)
    [i1,i2] = find(data.alternate(end).origin_time_loc==j,1,'first');
    if(data.alternate(end).origin_time_corr(i1,end)>=0)
        meta_origin2(j) = data.alternate(end).origin_time_corr(i1,end);
    else
        meta_origin2(j) = data.alternate(end).origin_time_corr(i1,i2);
    end
end
[~,mord2]=sort(meta_origin2);
orig_times2 = meta_origin2(mord2);
met_cumdist2 = cumsum(res.key.META_ORIGIN(mord2,end-1));

% calculate third compromise tumor count
min1 = min(orig_times1);
lmind1 = find(orig_times1==min1,1,'last');
min2 = min(orig_times2);
lmind2 = find(orig_times2==min2,1,'last');

if(min2>min1 || (min2==min1&&lmind2<lmind1))
    % append vestigial tumors from orig_times1 as extra tumors to
    % orig_times2 to make orig_times3
    inter = find(orig_times1==min2);
    if(isempty(inter))
        junc = find(orig_times1<min2,1,'last');
    else
        shift = max(0,length(inter)-lmind2);
        junc = inter(1)+shift-1;
    end
    orig_times = [orig_times1(1:junc); orig_times2];
    met_cumdist = [met_cumdist1(1:junc); met_cumdist1(junc)+met_cumdist2];
    mord = [mord1(1:junc); mord2];
elseif(min2<min1 || (min2==min1&&lmind2>lmind1))
    % append vestigial tumors from orig_times2 as extra tumors to
    % orig_times1 to make orig_times3
    inter = find(orig_times2==min1);
    if(isempty(inter))
        junc = find(orig_times2<min1,1,'last');
    else
        shift = max(0,length(inter)-lmind1);
        junc = inter(1)+shift-1;
    end
    orig_times = [orig_times2(1:junc); orig_times1];
    met_cumdist = [met_cumdist2(1:junc); met_cumdist2(junc)+met_cumdist1];
    mord = [mord2(1:junc); mord1];
else
    orig_times = orig_times2;
    met_cumdist = met_cumdist2;
    mord = mord2;
end


det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end);
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
% traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
% est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';
traj_full = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,res.t(end),res.capacity_factor);
ind_det_sizes = sum(traj_full<=det_sizes,2); % this finds the index of first time size becomes larger than det_size for each det_size; this is inefficient but don't want to write a loop
time_to_det_sizes = res.t3(ind_det_sizes)';

% calculate xcomp and difs based on birth times and size at detection
% to combine dif1 and dif2, we need to make sure the differences are in the
% same units/scales - otherwise one might be several magnitudes larger than
% the other. One solution is to convert the size difference represented by
% dif2 to a time difference using growth parameters
xcomp = met_cumdist;
dif1 = (orig_times-mod_births);
dif2 = (time_to_det_sizes-growth_times);
dif3 = [(orig_times(2:end)-orig_times(1:end-1))-(mod_births(2:end)-mod_births(1:end-1)); 0]; 

dif = sqrt(dif1.^2+dif2.^2+dif3.^2);

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [dif,xcomp] = getMetastasisDif_byCountAndSize3(res,data,i)
% Note: there's a situation that can occur do to the fact that we are only
% looking birth times and not tracking the size of metastasis at the time
% of detection. This means this fitter may give erroneous results
% 2/1/21: changing normalization of dif1 and dif2 since dividing by
% met_count makes it so that fitting low counts are empasized too much
% 3/2/21: accounting for number of metastasis at each data point

% calculate birth time step function
[~,mord]=sort(res.key.META_ORIGIN(:,end));
orig_times = res.key.META_ORIGIN(mord,end);
det_times = data.t(res.key.META_ORIGIN(mord,1))';
growth_times = det_times-orig_times;
det_sizes = res.key.META_ORIGIN(mord,3);
met_dist = res.key.META_ORIGIN(mord,end-1);
met_cumdist = cumsum(met_dist);
% calculate model tumor count
cum_meta = cumsum(res.num_meta);

% dif is defined as difference between origin times for each tumor so need
% to get origin times from the model based on times we get to integer tumor
% values
mod_births = zeros(length(orig_times),1); 
for j = 1:length(orig_times)
    locj = find(cum_meta>=met_cumdist(j),1);
    if(locj)
        mod_births(j) = res.t3(locj);
    else
        % set birth time to last time in simulation
        if(j>1)
            last_time = round(mod_births(j-1)/res.key.FIELD.dt)+1;
            time_per_met = res.key.FIELD.dt/res.num_meta(last_time);
        else
            time_per_met = inf;
        end
        if(j>1 && ~isinf(time_per_met))
            mod_births(j:end) = mod_births(j-1)+(met_cumdist(j:end)-met_cumdist(j-1))*time_per_met;
        else
            mod_births(j:end) = 100*res.t3(end);
        end
        break;
    end
end

% need to calculate expected size for each metastasis at detection time
% traj = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,ceil(max(growth_times)),res.capacity_factor);
% est_size_atDet = traj(round(growth_times/res.key.FIELD.dt)+1)';
traj_full = getSizeTrajectory_v3_4(res.key,1,res.key.FIELD.dt,res.t(end),res.capacity_factor);
ind_det_sizes = sum(traj_full<=det_sizes,2); % this finds the index of first time size becomes larger than det_size for each det_size; this is inefficient but don't want to write a loop
time_to_det_sizes = res.t3(ind_det_sizes)';

% calculate xcomp and difs based on birth times and size at detection
% to combine dif1 and dif2, we need to make sure the differences are in the
% same units/scales - otherwise one might be several magnitudes larger than
% the other. One solution is to convert the size difference represented by
% dif2 to a time difference using growth parameters
% -multiply by distribution magnitude to account for cases where multiple
% tumors are grouped together in the data
xcomp = met_cumdist;
dif1 = (orig_times-mod_births).*met_dist;
dif2 = (time_to_det_sizes-growth_times).*met_dist;
dif = sqrt(dif1.^2+dif2.^2);

% % make plot
% stepx = zeros(length(orig_times)*2,1);
% stepcdf = zeros(length(orig_times)*2,1);
% 
% stepx((1:length(orig_times))*2-1) = orig_times;
% stepx((1:length(orig_times))*2) = orig_times;
% stepcdf((1:length(orig_times))*2-1) = [0; met_cumdist(1:end-1)];
% stepcdf((1:length(orig_times))*2) = met_cumdist;
% 
% figure;
% plot(res.t3,cumsum(res.num_meta),stepx,stepcdf,'o-',abs(dif),xcomp,'--d')
% legend('Model','Data','Dif');
% xlabel('Time (Days)');
% ylabel('Number of Metastasis (Scaled)');
end

function [w] = getMMDWeights(res,data,i,xcomp,x_limit)
ixlim = res.key.FIELD.indexOf(x_limit,res.key);
dxs = res.key.FIELD.dxs;
w = zeros(1,length(dxs));
w(ixlim:end) = dxs(ixlim:end);
end

function [weight_func] = getMMDWeights_generator(x_limit)
    weight_func = @(res,data,i,xcomp) getMMDWeights(res,data,i,xcomp,x_limit);
end

function [w] = getPDFMMDWeights(res,data,i,xcomp,x_limit)
t2 = res.t2;
x = res.x';
dxs = res.key.FIELD.dxs';
datx = max(1,data.x(data.selector{i})');
datdist = data.dist(data.selector{i},i);

% construct approximate data distribution using Gaussian kernel
% g = @(x,xi,Di) (exp(-(x-xi).^2./Di/2)./sqrt(2*pi()*Di))*datdist;
cg = @(x,xi,Di) (1/2-erf((x-xi)./sqrt(2*Di))/2)*datdist;
gm = @(x,xi,Di) -diff(cg(x,xi,Di))./diff(x);

idatx = res.key.FIELD.indexOf(datx,res.key);
dx_datx = dxs(idatx)';
Di_corr = dx_datx.*2;
Di = Di_corr.*t2(i).*max(0.5*(res.key.RATES.growth(datx)+res.key.RATES.death(datx)+res.key.RATES.shed(datx)),1);

gmx1 = gm(x,datx,Di);
dist2 = res.dist2(:,i);
pdf = abs(gmx1)+abs(dist2);

% index of size limit
ixlim = res.key.FIELD.indexOf(x_limit,res.key);

% weights are proportional to estimated data distribution
w = zeros(1,length(dxs));
w(ixlim:end) = pdf(ixlim:end).*dxs(ixlim:end)/sum(pdf(ixlim:end))*length(dxs(ixlim:end));
end

function [weight_func] = getPDFMMDWeights_generator(x_limit)
    weight_func = @(res,data,i,xcomp) getPDFMMDWeights(res,data,i,xcomp,x_limit);
end

function [w] = noWeights(res,data,i,xcomp)
w = ones(1,length(xcomp))/length(xcomp);
end

function [w] = distanceWeights(res,data,i,xcomp)
da = log10(max(data.a(data.selector{i}),1));
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

function [w] = pdfWeights(res,data,i,xcomp)
dx = data.x(data.selector{i});
dcdf = data.cdf(data.selector{i},i);
del = dx(2:end)-dx(1:end-1);
dis_ave = (dcdf(2:end)-dcdf(1:end-1))./del;
w = 1./dis_ave.^2;
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

function [w] = noWeights_TimeDependentLowFilter(res,data,i,xcomp,limits)
da = data.a(data.selector{i});
da_mid = 0.5*(da(2:end)+da(1:end-1));
w = ones(size(da_mid));
w(da(1:end-1)<limits(i)) = 0;
if(~isempty(w))
    w(end) = 1; % make sure at least one nonzero val
end
w = w/sum(w);
if(iscolumn(w))
    w = w';
end

end

function [weight_func] = noWeights_TimeDependentLowFilter_generator(start_limit)
    traj_func = @(res,data,x0) getSizeTrajectory_v3_4(res.key,x0,res.key.FIELD.dt,data.t(end)-data.t(1),res.capacity_factor);
    subsel = @(vec,sel) vec(sel);
    weight_func = @(res,data,i,xcomp) noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
    subsel(traj_func(res,data,start_limit),round((data.t-data.t(1))/res.key.FIELD.dt)+1)));
end

function [w] = noWeights_PrimaryFilter(res,data,i,xcomp)
% find size of smallest primary (when ccdf==ccdf(1,1) or if
% ccdf<ccdf(1,1), the smallest tumor size in the system)
w = noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
   getPrimaryEndSizes(data)));
end

function [w] = noWeights_MetaFilter(res,data,i,xcomp)
% find size of smallest primary (when ccdf==ccdf(1,1) or if
% ccdf<ccdf(1,1), the smallest tumor size in the system)
% in this case, the weights for primaries are subtracted so only metastasis
% contribute to weights
w = max(0,noWeights(res,data,i,xcomp)-noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
   getPrimaryEndSizes(data))));
if(sum(w)==0)
    w = ones(size(w));
end
w = w/sum(w);
end

function [p_ends] = getPrimaryEndSizes(data)
% version using data
p_ends = zeros(1,length(data.t));
dcdf = data.cdf(data.selector{1},1);
if(isfield(data,'cdf_adj'))
    p_adj = data.cdf_adj(1);
else
    p_adj = 0;
end
dccdf1 = p_adj+dcdf(end);
for i = 1:length(data.t)
    dcdf = data.cdf(data.selector{i},i);
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(i);
    else
        p_adj = 0;
    end
    dccdfi = p_adj+dcdf(end)-[0; dcdf(1:end-1)];
    p_ei = find(dccdfi<=dccdf1,1,'first');
    p_ends(i) = data.x(data.selector{i}(p_ei));
end
% % version using model results
% ccdf_p = (res.dist_cum2(res.max_index+1,:)-res.dist_cum2(1:res.max_index+1,:))<=(res.dist_cum2(res.max_index+1,1)-res.dist_cum2(1,1));
% p_inds = zeros(1,size(ccdf_p,2));
% for i = 1:size(ccdf_p,2)
%     p_inds(i) = max(find(ccdf_p(:,i),1,'first')-1,1);
% end
% p_ends = res.x(p_inds);
end

function [w] = heightAdjustedWeights(res,data,i,xcomp)
dcdf = data.cdf(data.selector{i},i);
if(isfield(data,'cdf_adj'))
    p_adj = data.cdf_adj(i);
else
    p_adj = 0;
end
dccdf1 = p_adj+dcdf(end)-[0; dcdf(1:end-1)];

avg_height = 0.5*(dccdf1(2:end)+dccdf1(1:end-1));
denom = avg_height;
w = 1./denom.^2;
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end


function [weight_func] = heightAdjustedWeights_lowFilter_generator(limit)
    weight_func1 = @(res,data,i,xcomp) heightAdjustedWeights(res,data,i,xcomp).*...
        noWeights_lowFilter(res,data,i,xcomp,res.key.CONV.v2a(limit));
    weight_func = @(res,data,i,xcomp) weight_func1(res,data,i,xcomp)/sum(weight_func1(res,data,i,xcomp));
end

function [weight_func] = heightAdjustedWeights_TimeDependentLowFilter_generator(start_limit)
    traj_func = @(res,data,x0) getSizeTrajectory_v3_4(res.key,x0,res.key.FIELD.dt,data.t(end)-data.t(1),res.capacity_factor);
    subsel = @(vec,sel) vec(sel);
    weight_func1 = @(res,data,i,xcomp) heightAdjustedWeights(res,data,i,xcomp).*...
        noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
    subsel(traj_func(res,data,start_limit),round((data.t-data.t(1))/res.key.FIELD.dt)+1)));
    weight_func = @(res,data,i,xcomp) weight_func1(res,data,i,xcomp)/sum(weight_func1(res,data,i,xcomp));
end

function [weight_func] = distanceWeights_TimeDependentLowFilter_generator(start_limit)
    traj_func = @(res,data,x0) getSizeTrajectory_v3_4(res.key,x0,res.key.FIELD.dt,data.t(end)-data.t(1),res.capacity_factor);
    subsel = @(vec,sel) vec(sel);
    weight_func1 = @(res,data,i,xcomp) distanceWeights(res,data,i,xcomp).*...
        noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
    subsel(traj_func(res,data,start_limit),round((data.t-data.t(1))/res.key.FIELD.dt)+1)));
    weight_func = @(res,data,i,xcomp) weight_func1(res,data,i,xcomp)/sum(weight_func1(res,data,i,xcomp));
end

function [weight_func] = distanceWeights_lowFilter_generator(limit)
    weight_func1 = @(res,data,i,xcomp) distanceWeights(res,data,i,xcomp).*...
        noWeights_lowFilter(res,data,i,xcomp,res.key.CONV.v2a(limit));
    weight_func = @(res,data,i,xcomp) weight_func1(res,data,i,xcomp)/sum(weight_func1(res,data,i,xcomp));
end

function [weight_func] = pdfWeights_lowFilter_generator(limit)
    weight_func1 = @(res,data,i,xcomp) pdfWeights(res,data,i,xcomp).*...
        noWeights_lowFilter(res,data,i,xcomp,res.key.CONV.v2a(limit));
    weight_func = @(res,data,i,xcomp) weight_func1(res,data,i,xcomp)/sum(weight_func1(res,data,i,xcomp))*data.cdf(data.selector{i}(end),i);
end

function [w] = distanceWeights_PrimaryFilter(res,data,i,xcomp)
% find size of smallest primary (when ccdf==ccdf(1,1) or if
% ccdf<ccdf(1,1), the smallest tumor size in the system)
w = distanceWeights(res,data,i,xcomp).*...
    noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
   getPrimaryEndSizes(data)));
end

function [w] = pdfWeights_PrimaryFilter(res,data,i,xcomp)
% find size of smallest primary (when ccdf==ccdf(1,1) or if
% ccdf<ccdf(1,1), the smallest tumor size in the system)
w = pdfWeights(res,data,i,xcomp).*...
    noWeights_TimeDependentLowFilter(res,data,i,xcomp,res.key.CONV.v2a(...
   getPrimaryEndSizes(data)));
w = w/sum(w)*data.cdf(data.selector{i}(end),i);
end

function [w] = countWeights(res,data,i,xcomp)
dcdf = data.cdf(data.selector{i},i);
w = dcdf(2:end)-dcdf(1:end-1);
w = w/sum(w);
if(iscolumn(w))
    w = w';
end
end

function [plotter] = getPlotter()
    plotter = struct();
    plotter.start = @(res,data,s) figure;
    plotter.continue = @plotterContinue3;
    plotter.chain = @plotterChain;
    plotter.end = @plotterEnd;
    plotter.chainEnd = @plotterChainEnd;
end

function [plotter] = getPlotterNoInference()
    plotter = struct();
    plotter.start = @(res,data,s) figure;
    plotter.continue = @plotterContinue2;
    plotter.chain = @plotterChain;
    plotter.end = @plotterEndNoInference;
    plotter.chainEnd = @plotterChainEnd;
end

function [] = plotterContinue(plot_func,res,data,i,dif,xcomp)
if(any(i==data.plot_indices))
    da = data.x(data.selector{i});
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(i);
    else
        p_adj = 0;
    end
    if(~isempty(da))
        % dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
        dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
        imin = 1;
        imax = res.max_index;
        ma = res.x(imin:imax+1);
        mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
        colors = hsv(length(data.t));
        dac = zeros(2*length(da),1);
        dac(1:2:end-1) = da;
        dac(2:2:end) = da;
        dccdfc = zeros(2*length(da),1);
        dccdfc(1:2:end-1) = dccdf1;
        dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
        plot_func(dac,dccdfc,'o:',ma,mccdf,'-',xcomp,abs(dif)/max(abs(dif)),':',...
            'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
            'LineWidth',1.5,'MarkerSize',4);
        hold on;
    else
        plot_func(1,0,1,0,1,0);
        hold on;
    end
end
end

function [] = plotterContinue2(plot_func,res,data,difs,xcomps)
for i = 1:length(data.t)
    dif = difs{i};
    xcomp = xcomps{i};
    if(any(i==data.plot_indices))
        da = data.x(data.selector{i});
        if(isfield(data,'cdf_adj'))
            p_adj = data.cdf_adj(i);
        else
            p_adj = 0;
        end
        
        if(~isempty(da))
            % dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
            dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
            imin = 1;
            imax = res.max_index;
            ma = res.x(imin:imax+1);
            mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
            colors = hsv(length(data.t));
            dac = zeros(2*length(da),1);
            dac(1:2:end-1) = da;
            dac(2:2:end) = da;
            dccdfc = zeros(2*length(da),1);
            dccdfc(1:2:end-1) = dccdf1;
            dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
            plot_func(dac,dccdfc,'o:',ma,mccdf,'-',xcomp,abs(dif)/max(abs(dif)),'-.',...
                'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
                'LineWidth',1.5,'MarkerSize',4);
            hold on;
        else
            plot_func(1,0,1,0,1,0);
            hold on;
        end
    end
end
end

function [] = plotterContinue3(plot_func,res,data,difs,xcomps)
for i = 1:length(data.t)
    dif = difs{i};
    xcomp = xcomps{i};
    if(any(i==data.plot_indices))
        da = data.x(data.selector{i});
        if(isfield(data,'cdf_adj'))
            p_adj = data.cdf_adj(i);
        else
            p_adj = 0;
        end
        
        if(~isempty(da))
            % dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
            %             dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
            imin = 1;
            imax = res.max_index;
            ma = res.x(imin:imax+1);
            mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
            colors = hsv(length(data.t));
            %             dac = zeros(2*length(da),1);
            %             dac(1:2:end-1) = da;
            %             dac(2:2:end) = da;
            %             dccdfc = zeros(2*length(da),1);
            %             dccdfc(1:2:end-1) = dccdf1;
            %             dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
            dac = data.alternate(i).step_ccdf(:,1);
            dccdfc = data.alternate(i).step_ccdf(:,2);
            dac1 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,1);
            dccdfc1 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,2);
            if(isempty(dac1))
                dac1 = 1;
                dccdfc1 = 0;
            end
            dac2 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,1);
            dccdfc2 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,2);
            if(isempty(dac2))
                dac2 = 1;
                dccdfc2 = 0;
            end
            %             plot_func(dac,dccdfc,'--',ma,mccdf,'-',dac1,dccdfc1,'s',...
            %                 dac2,dccdfc2,'d',xcomp,abs(dif)/max(abs(dif)),'+',...
            %                 'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
            %                 'LineWidth',1.5,'MarkerSize',4);
            plot_func(dac,dccdfc,'--',ma,mccdf,'-',dac1,dccdfc1,'s',...
                dac2,dccdfc2,'d',...
                'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
                'LineWidth',1.5,'MarkerSize',4);
            hold on;
        else
            plot_func(1,0,1,0,1,0);
            hold on;
        end
    end
end
end

function [] = plotterChain(plot_func,res,data,chains,xcomps)
for i = 1:length(data.t)
    xcomp = xcomps{i};
    if(any(i==data.plot_indices))
        da = data.x(data.selector{i});
        if(isfield(data,'cdf_adj'))
            p_adj = data.cdf_adj(i);
        else
            p_adj = 0;
        end
%         dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
        dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
        imin = 1;
        imax = res.max_index;
        ma = res.x(imin:imax+1);
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
            chain = chains{j,i};
            if(any(chain))
                ylim = max([data.cdf(end,:) res.dist_cum2(end,:)]);
                ylim = 10*floor(0.1*ylim)+10;
                chain = chain/max(abs(chain))*ylim/2+ylim/2;
                plot_func(xcomp,chain,'--',...
                    'Color',max(0,(0.6-0.05*j))*colors(i,:));
            end
        end
    end
end
end

function [] = plotterEnd(res,data,fit)
hold off;
title(res.key.TITLE,'fontsize',7);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
% get legend
leg = cell(1,4*length(data.plot_indices));
for i = 1:length(data.plot_indices)
    leg{4*i-3} = ['Data, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{4*i-2} = ['Model, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{4*i-1} = ['Inf. Prim, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{4*i-0} = ['Inf. Meta, t = ' num2str(data.t(data.plot_indices(i)))];
    %     leg{5*i} = ['fit = ' sprintf('%1.3f',fit(data.plot_indices(i)))];
end
legend(leg,'Location','Best');
alow = data.lowa;
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.x > xlow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
% ylim = max([data.cdf(end,:)-data.cdf(ia1,:) res.dist_cum2(end,:)-res.dist_cum2(ix1,:)]);
% ylim = floor(ylim)+1; 
ylim = max([data.cdf(end,:) res.dist_cum2(end,:)]);
ylim = 10*floor(0.1*ylim)+10;
axis([alow inf 0.01 ylim]);
end

function [] = plotterEndNoInference(res,data,fit)
hold off;
title(res.key.TITLE,'fontsize',7);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
% get legend
leg = cell(1,2*length(data.plot_indices));
for i = 1:length(data.plot_indices)
    leg{2*i-1} = ['Data, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{2*i} = ['Model, t = ' num2str(data.t(data.plot_indices(i)))];
end
legend(leg,'Location','Best');
alow = data.lowa;
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.x > xlow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
% ylim = max([data.cdf(end,:)-data.cdf(ia1,:) res.dist_cum2(end,:)-res.dist_cum2(ix1,:)]);
% ylim = floor(ylim)+1; 
ylim = max([data.cdf(end,:) res.dist_cum2(end,:)]);
ylim = 10*floor(0.1*ylim)+10;
axis([alow inf 0.01 ylim]);
end

function [] = plotterChainEnd(res,data,grad,chains)
hold off;
title(res.key.TITLE,'fontsize',7);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
nplot = length(data.plot_indices);
% get legend
% find which chains have nonzero values
nl0 = 2*nplot;
for i = 1:size(chains,2)
    for j = 1:size(chains,1)
        nl0 = nl0+any(chains{j,i});
    end
end
nl = length(chains)+2;
parnames = fieldnames(res.key.PARAMETERS);
leg = cell(1,nl0);
nlc = 0;
for i = 1:length(data.plot_indices)
    leg{nlc+1} = ['Data, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{nlc+2} = ['Model, t = ' num2str(data.t(data.plot_indices(i)))];
    nlc = nlc+2;
    for j = 3:nl
        if(any(chains{j-2,i}))
            nlc = nlc+1;
            leg{nlc} = [parnames{j-2} ': <' num2str(mean(chains{j-2,i}))...
                ', ' num2str(grad(j-2,i)) '>'];
        end
    end
end
legend(leg,'Location','Best');
alow = 1;
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.x > xlow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
% ylim = max([data.cdf(end,:)-data.cdf(ia1,:) res.dist_cum2(end,:)-res.dist_cum2(ix1,:)]);
% ylim = floor(ylim)+1;
ylim = max([data.cdf(end,:) res.dist_cum2(end,:)]);
ylim = 10*floor(0.1*ylim)+10;
axis([alow inf 0 ylim]);
end

function [plotter] = getLinkedPlotter()
    plotter = struct();
    plotter.start = @linkedPlotterStart;
    plotter.continue = @linkedPlotterContinue2;
    plotter.end = @linkedPlotterEnd;
    plotter.chain = @linkedPlotterContinue2;
    plotter.chainEnd = @(res,data,grad,chains) linkedPlotterEnd(res,data,grad);
end

function [] = linkedPlotterStart(res,data,s)
if(evalin('base','exist(''fitter_figure'',''var'')'))
    if(~evalin('base','fitter_figure.isvalid'))
        evalin('base','fitter_figure = figure;');
        initializeLinkedPlot2(s.plot,res,data);
    end
else
    evalin('base','fitter_figure = figure;');
    initializeLinkedPlot2(s.plot,res,data);
end
end

function [] = initializeLinkedPlot(plot_func,res,data)
% initialize plot
evalin('base','fitter_figure.Position = [1 1 950 680];');
% update 5/19/20: plotting volume
assignin('base','fitter_rx',res.x(1:end));
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
    dat_ai = data.x(data.selector{i});
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
    '(fitter_rx,fitter_res,''LineWidth'',1.25);']);
hold on;
for i = 1:length(data.t)
    evalin('base',[func2str(plot_func) '(fitter_da{' num2str(i) '},'...
        'fitter_dat{' num2str(i) '},''o--'',''LineWidth'',0.75);']);
end
hold off;
title(res.key.TITLE,'fontsize',7);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
evalin('base','linkdata(fitter_figure)');

% store fit variable
evalin('base','fitter_min_fit = inf;');
end

function [] = initializeLinkedPlot2(plot_func,res,data)
% uses data.alternate structure (getKey_backtrack_v5_1 and later)
% initialize plot
evalin('base','fitter_figure.Position = [1 1 950 680];');
% update 5/19/20: plotting volume
assignin('base','fitter_rx',res.x(1:end));
assignin('base','fitter_res',addByColumn(...
    -res.dist_cum2(1:end,:),res.dist_cum2(end,:)));

for i = 1:length(data.t)
    assignin('base',['fitter_da' num2str(i)],data.alternate(i).step_ccdf(:,1));
    assignin('base',['fitter_dat' num2str(i)],data.alternate(i).step_ccdf(:,2));
    assignin('base',['fitter_ipx' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,1));
    assignin('base',['fitter_ipcdf' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,2));
    assignin('base',['fitter_imx' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,1));
    assignin('base',['fitter_imcdf' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,2));
end

% make plots
colors = hsv(length(data.t));
plotcol = [colors;0.6*colors];
assignin('base','fitter_colors',plotcol);
evalin('base','figure(fitter_figure);');
evalin('base','plot(0,0);');
evalin('base','set(fitter_figure.CurrentAxes,''ColorOrder'',fitter_colors);');
evalin('base','set(fitter_figure.CurrentAxes,''NextPlot'',''replacechildren'');');

data_string1 = '';
for i = 1:length(data.t)-1
    data_string1 = [data_string1 'fitter_da' num2str(i) ','...
        'fitter_dat' num2str(i) ',''--'','];
end
i = length(data.t);
data_string1 = [data_string1 'fitter_da' num2str(i) ','...
    'fitter_dat' num2str(i) ',''--'''];

data_string2 = '';
for i = 1:length(data.t)-1
    data_string2 = [data_string2 'fitter_ipx' num2str(i) ','...
        'fitter_ipcdf' num2str(i) ',''s'',fitter_imx' num2str(i) ','...
        'fitter_imcdf' num2str(i) ',''d'','];
end
data_string2 = [data_string2 'fitter_ipx' num2str(i) ','...
        'fitter_ipcdf' num2str(i) ',''s'',fitter_imx' num2str(i) ','...
        'fitter_imcdf' num2str(i) ',''d'''];
    
evalin('base',[func2str(plot_func) ...
    '(fitter_rx,fitter_res,' data_string1 ',' data_string2 ',''LineWidth'',1.25);']);
    
title(res.key.TITLE,'fontsize',7);
xlabel(data.xlabel);
ylabel('Cumulative Number of Tumors');
evalin('base','linkdata(fitter_figure)');

% store fit variable
evalin('base','fitter_min_fit = inf;');
end

function [] = linkedPlotterContinue(plot_func,res,data,difs,xcomps)
% 5/19/20: plotting with volume...
assignin('base','fitter_rx',res.x(1:end));
assignin('base','fitter_res',addByColumn(...
    -res.dist_cum2(1:end,:),res.dist_cum2(end,:)));
end

function [] = linkedPlotterContinue2(plot_func,res,data,difs,xcomps)
% 5/19/20: plotting with volume...
assignin('base','fitter_rx',res.x(1:end));
assignin('base','fitter_res',addByColumn(...
    -res.dist_cum2(1:end,:),res.dist_cum2(end,:)));

for i = 1:length(data.t)
    assignin('base',['fitter_da' num2str(i)],data.alternate(i).step_ccdf(:,1));
    assignin('base',['fitter_dat' num2str(i)],data.alternate(i).step_ccdf(:,2));
    assignin('base',['fitter_ipx' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,1));
    assignin('base',['fitter_ipcdf' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,2));
    assignin('base',['fitter_imx' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,1));
    assignin('base',['fitter_imcdf' num2str(i)],data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,2));
end
end

function [] = linkedPlotterEnd(res,data,fit)
evalin('base','figure(fitter_figure);');
title(res.key.TITLE,'fontsize',7);
% get legend
leg = cell(1,2*length(data.t));
for i = 1:length(data.t)
    leg{i+length(data.t)} = ['Data, t = ' num2str(data.t(i))];
    leg{i} = ['Model Fit =' sprintf(' %1.3f',fit(:,i))];
end
legend(leg,'Location','Best');
alow = 10^floor(log10(max(data.lowx,data.x(1))));
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.x > alow,1);
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
axis([alow inf 0.01 ylim]);
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

function [plotter] = getDifPlotter()
    plotter = struct();
    plotter.start = @(res,data,s) figure;
    plotter.continue = @difPlotterContinue;
    plotter.chain = @plotterChain;
    plotter.end = @difPlotterEnd;
    plotter.chainEnd = @plotterChainEnd;
end

function [] = difPlotterContinue(plot_func,res,data,difs,xcomps)
for i = 1:length(data.t)
    dif = difs{i};
    xcomp = xcomps{i};
    if(any(i==data.plot_indices))
        da = data.x(data.selector{i});
        if(isfield(data,'cdf_adj'))
            p_adj = data.cdf_adj(i);
        else
            p_adj = 0;
        end
        
        if(~isempty(da))
            % dccdf = p_adj+data.cdf(data.selector{i}(end),i)-data.cdf(data.selector{i},i);
            %             dccdf1 = p_adj+data.cdf(data.selector{i}(end),i)-[0; data.cdf(data.selector{i}(1:end-1),i)];
            imin = 1;
            imax = res.max_index;
            ma = res.x(imin:imax+1);
            mccdf = res.dist_cum2(imax+1,i)-res.dist_cum2(imin:imax+1,i);
            colors = hsv(length(data.t));
            %             dac = zeros(2*length(da),1);
            %             dac(1:2:end-1) = da;
            %             dac(2:2:end) = da;
            %             dccdfc = zeros(2*length(da),1);
            %             dccdfc(1:2:end-1) = dccdf1;
            %             dccdfc(2:2:end) = [dccdf1(2:end); p_adj];
            dac = data.alternate(i).step_ccdf(:,1);
            dccdfc = data.alternate(i).step_ccdf(:,2);
            dac1 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,1);
            dccdfc1 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==1,2);
            if(isempty(dac1))
                dac1 = 1;
                dccdfc1 = 0;
            end
            dac2 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,1);
            dccdfc2 = data.alternate(i).step_ccdf(data.alternate(i).step_ccdf(:,3)==2,2);
            if(isempty(dac2))
                dac2 = 1;
                dccdfc2 = 0;
            end
            %             plot_func(dac,dccdfc,'--',ma,mccdf,'-',dac1,dccdfc1,'s',...
            %                 dac2,dccdfc2,'d',xcomp,abs(dif)/max(abs(dif)),'+',...
            %                 'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
            %                 'LineWidth',1.5,'MarkerSize',4);
            %             plot_func(dac,dccdfc,'--',ma,mccdf,'-',dac1,dccdfc1,'s',...
            %                 dac2,dccdfc2,'d',...
            %                 'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
            %                 'LineWidth',1.5,'MarkerSize',4);
            %             plot_func(xcomp,dif,'-o','Color',0.6*colors(i,:),'LineWidth',1.25)
            plot_func(dac,dccdfc,'--',ma,mccdf,'-',...
                xcomp,abs(dif),'d:',...
                'MarkerEdgeColor',colors(i,:),'Color',0.6*colors(i,:),...
                'LineWidth',1.5,'MarkerSize',4);
            hold on;
        else
            plot_func(1,0,1,0,1,0);
            hold on;
        end
    end
end
end

function [] = difPlotterEnd(res,data,fit)
hold off;
title(res.key.TITLE,'fontsize',7);
xlabel(data.xlabel);
ylabel('Weighted Difference');
% get legend
leg = cell(1,3*length(data.plot_indices));
for i = 1:length(data.plot_indices)
    leg{3*i-2} = ['Dat, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{3*i-1} = ['Mod, t = ' num2str(data.t(data.plot_indices(i)))];
    leg{3*i} = ['Dif, t = ' num2str(data.t(data.plot_indices(i)))];
end
legend(leg,'Location','Best');
alow = data.lowa;
xlow = res.key.CONV.a2v(alow);
ia1 = find(data.x > xlow,1);
if(isempty(ia1))
    ia1 = 1;
end
ix1 = find(res.x > xlow,1);
if(isempty(ix1))
    ix1 = 1;
end
% ylim = max([data.cdf(end,:)-data.cdf(ia1,:) res.dist_cum2(end,:)-res.dist_cum2(ix1,:)]);
% ylim = floor(ylim)+1; 
% ylim = max([data.cdf(end,:) res.dist_cum2(end,:)]);
% ylim = 10*floor(0.1*ylim)+10;
axis([alow inf -inf inf]);
end

function [chains] = getFullChain(res,data,xcomps)
chains = getFullChain_v2_0(res,data,xcomps);
end

function [chains] = getBlankChain(res,data,xcomps)
key = res.key;
pars = key.PARAMETERS;
parnames = fieldnames(key.PARAMETERS);
npar = length(parnames);
rates = key.RATES;
f0 = rates.size1_ratio;
nt = length(data.t);
chains = cell(npar,nt);

for i = 1:npar
    for j = 1:nt
        chains{i,j} = zeros(length(xcomps{j}),1);
    end
end
end