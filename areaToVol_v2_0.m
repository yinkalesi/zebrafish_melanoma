%% areaToVol_v2_0
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 3/24/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%  areaToVol is meant to convert cell area to volume
%  area: float; cell area from experiment
%  vol: float; approximate cell volume
%  [vol] = areaToVol_v2_0(area)
%% Version History
%  2.0: separating fL_cell into multiple constants to clarify and adding a
%  thickness limit to tumors are assumed to be more of an ellipsoid

function [vol] = areaToVol_v2_0(area)

um_px = (4.22/1388+3.16/1040)/2*1e4; %alternatively, um_px = sqrt(4.22*3.16*(1e8)/1388/1040)
fL_cell =  4/3*pi*5^3; % diameter of cancer cell = 10 um
shape_factor = 4/3/sqrt(pi); % value assumes sphere
calib_factor = 1/10; % *needs calibration
cell_px3 = um_px^3/fL_cell;
thick_limit = 0.4*1e4/um_px; % assuming tumor thickness limit of 0.5 cm

thick_adj = 1./(1+sqrt(area/pi)/thick_limit);
vol = cell_px3*shape_factor*calib_factor*(area).^(3/2).*thick_adj;
