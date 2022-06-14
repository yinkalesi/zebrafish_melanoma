%% getPixelLine_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/19/19
%  Project: Fish Image Analysis
%  get line of pixels between two pixels

function [line] = getPixelLine_v1_0(p1,p2)

x = p2(1)-p1(1);
y = p2(2)-p1(2);
n = max(ceil(abs(x)),ceil(abs(y)));
dx = x/n;
dy = y/n;

px = [p1(1) floor(p1(1)+(1:n-1)*dx) p2(1)];
py = [p1(2) floor(p1(2)+(1:n-1)*dy) p2(2)]; 

line = [px' py'];