%% getConverter
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 3/27/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%  converter: structure with conversion functions
%  [converter] = getConverter()
%% Version History
%  2.0: update to account for tumor thickness, calibrationn

function [converter] = getConverter()

converter = struct();
converter.a2v = @areaToVol_v2_0;
converter.v2a = @volToArea_v2_0;