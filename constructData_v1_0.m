function [edata] = constructData_v1_0(t,sizes,dist,converter)
%% constructData_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 6/17/19
%  Project: Tumor Growth, Logarithmic Continuum Form
%  specifically for zebrafish data

% select which data to use
%% Version History
%  1.0: minor derived from getExperimentalData_v3_0; want to be able to
%  create arbitrary data structs for testing purposes
%  3/24/21: for size=0, set dist to 0

%% LOAD EXPERIMENTAL DATA

edata.dist = dist;
edata.dist(sizes==0,:) = 0; % assuming size = 0 is a placeholder
edata.a = sizes;
edata.t = t;
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

edata.conv = converter;

% determine data labels here and other plot options
edata.xlabel = 'Tumor Size';
edata.ylabel = 'Number of Tumors';

edata.plot_indices = 1:length(edata.t); % not used for every plot
% parameter for fitting
edata.weights = ones(1,length(edata.t));