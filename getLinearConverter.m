%% getLinearConverter
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 12/5/16
%  Project: Tumor Growth, Logarithmic Continuum Form
%  converter: structure with conversion functions
%  this converter does nothing
%  [converter] = getConverter()

function [converter] = getLinearConverter(dvda)

converter = struct();
% from data: dvda ~ 90.656511309673434
converter.a2v = @(x) dvda*x;
converter.v2a = @(x) x/dvda;