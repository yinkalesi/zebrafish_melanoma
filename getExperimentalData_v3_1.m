function [edata] = getExperimentalData_v3_1(fname,converter,picks)
%% getExperimentalData
%  Version 3.1
%  Author: Adeyinka Lesi
%  Date: 7/30/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%  specifically for zebrafish data

% select which data to use
%% Version History
%  3.0: minor changes for using intensity data instead
%  3.1: time shift (one is zero 3/30/21)

%% LOAD EXPERIMENTAL DATA
temp = struct();
[temp.dist, temp.a, temp.t] = loadDistribution(fname);
if(~exist('picks','var'))
    picks = 1:length(temp.t);
end
edata.dist = temp.dist(:,picks);
edata.a = temp.a;
edata.t = temp.t(picks)-1;
edata.x = converter.a2v(edata.a); % estimate of tumor size in cells
edata.cdf = zeros(length(edata.a),length(edata.t));
edata.selector = cell(1,length(edata.t));
for i = 1:length(edata.t)
    edata.cdf(:,i) = getStepIntegral(edata.a,edata.dist(:,i));
    edata.selector{i} = find(edata.dist(:,i)>0);    
end
edata.max = max(edata.x);

edata.lowa = converter.v2a(1); % lowest tumor area that can be reliably measured  
edata.lowx = converter.a2v(edata.lowa); % lowest cell size that can be reliably measured 

% % data below the threshold comes in with a size value of zero - model can't
% % handle that so we will shift it to center of region below lowa
% if(edata.a(1) == 0)
%     edata.x(1) = 1;
%     edata.a(1) = converter.v2a(edata.x(1));
% end
edata.conv = converter;

% determine data labels here and other plot options
edata.xlabel = 'Tumor Size';
edata.ylabel = 'Number of Tumors per Body';

edata.plot_indices = 1:length(edata.t); % not used for every plot
% parameter for fitting
edata.weights = ones(1,length(edata.t));