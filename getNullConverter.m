%% getNullConverter
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 3/23/16
%  Project: Tumor Growth, Logarithmic Continuum Form
%  converter: structure with conversion functions
%  this converter does nothing
%  [converter] = getConverter()

function [converter] = getNullConverter()

converter = struct();
converter.a2v = @(x) x;
converter.v2a = @(x) x;