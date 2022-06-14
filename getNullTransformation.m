%% getNullTransformation
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 4/17/16
%  Project: Tumor Growth, Logarithmic Continuum Form
%  getTransformationRule returns a set of functions to guide the spacing of
%  the mesh in the simulation
%  transform: struct; contains transform, inverse transform, transform
%  derivative and double derivative

function [transform] = getNullTransformation()

transform = struct();

transform.x2y = @(x) x;
transform.y2x = @(y) y;
transform.dydx = @(x) 1;
transform.d2ydx2 = @(x) 0;


