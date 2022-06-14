function [edata] = getExperimentalData_v3_2(fname,converter,picks)
%% getExperimentalData
%  Version 3.0
%  Author: Adeyinka Lesi
%  Date: 3/16/21
%  Project: Tumor Growth, Logarithmic Continuum Form
%  specifically for zebrafish data

% select which data to use
%% Version History
%  3.0: minor changes for using intensity data instead
%  3.1: time shift
%  3.2: size zero tumors, switched to FirstIsZero mode

%% LOAD EXPERIMENTAL DATA
temp = struct();
[temp.dist, temp.a, temp.t] = loadDistribution(fname);
if(~exist('picks','var'))
    picks = 1:length(temp.t);
end
edata.dist = temp.dist(:,picks);
edata.a = temp.a;
edata.t = temp.t(picks)-temp.t(1);
edata.x = converter.a2v(edata.a); % estimate of tumor size in cells

%    -0.2372    1.8632
%    -0.2191    2.3888
% replace tumors with size zero with extrapolated distribution
% extrapolation results day1: slope=-0.2372, ntums=1.8632
% linkage point -> (log10(xf) = 2, y = 1.394)
% data ntums fin -> y = 1.422
% tumor increment = 1/109;
% extrapolation results day2: slope=-0.2191, ntums=2.3888
% linkage point -> (log10(xf) = 2.2, y = 1.905)
% data ntums fin -> y = 1.945
% tumor increment = 1/109;

% day1 extrapolation
% d1_slope = -0.2372;
d1_add = edata.dist(1,1);
tum_inc1 = min(nonzeros(edata.dist(:,1)));%1/109;
% ntum_add1 = round(d1_add/tum_inc1);
% delx1 = -tum_inc1/d1_slope;
% x_add1 = 10.^((0:ntum_add1-1)*delx1)';
% a_add1 = converter.v2a(x_add1);
% dist_add1 = ones(size(x_add1))*tum_inc1;

[x_add1,dist_add1] = getExtrapolatedDistribution_Exp_v1_0(d1_add,tum_inc1,[1 100]);
a_add1 = converter.v2a(x_add1);

% day2 extrapolation
% d2_slope = -0.2191;
d2_add = 0;%40/109; %edata.dist(1,2);
tum_inc2 = min(nonzeros(edata.dist(:,2)));%1/109;
% ntum_add2 = round(d2_add/tum_inc2);
% delx2 = -tum_inc2/d2_slope;
% x_add2 = 10.^((0:ntum_add2-1)*delx2)';
% a_add2 = converter.v2a(x_add2);
% dist_add2 = ones(size(x_add2))*tum_inc2;

[x_add2,dist_add2] = getExtrapolatedDistribution_Exp_v1_0(d2_add,tum_inc2,[1 100]);
a_add2 = converter.v2a(x_add2);

new_a = [edata.a; a_add1; a_add2];
new_x = [edata.x; x_add1; x_add2];
new_dist = [edata.dist; [dist_add1 zeros(length(dist_add1),2)]; [zeros(length(dist_add2),1) dist_add2 zeros(length(dist_add2),1)]];
extrapol_label = [zeros(size(edata.x)); ones(size(x_add1)); 2*ones(size(x_add2))];
if(new_a(1)==0) % remove density due to size zero tumors on first and second days
    new_dist(1,1) = 0;
    new_dist(1,2) = 0;
end

[~,ord] = sort(new_a);
edata.a = new_a(ord);
edata.x = new_x(ord);
edata.dist = new_dist(ord,:);
edata.extrapol_label = extrapol_label(ord,:);

edata.cdf = zeros(length(edata.a),length(edata.t));
edata.selector = cell(1,length(edata.t));
for i = 1:length(edata.t)
    edata.cdf(:,i) = getStepIntegral(edata.a,edata.dist(:,i));
    edata.selector{i} = find(edata.dist(:,i)>0);    
end
edata.max = max(edata.x);

edata.lowa = converter.v2a(1); % lowest tumor area that can be reliably measured  
edata.lowx = converter.a2v(edata.lowa); % lowest cell size that can be reliably measured 

edata.conv = converter;

% determine data labels here and other plot options
edata.xlabel = 'Tumor Size';
edata.ylabel = 'Number of Tumors per Body';

edata.plot_indices = 1:length(edata.t); % not used for every plot
% parameter for fitting
edata.weights = ones(1,length(edata.t));