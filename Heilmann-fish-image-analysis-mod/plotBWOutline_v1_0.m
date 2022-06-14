%% plotBWOutline_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 02/04/20
%  Project: Image analysis
%% Version History
%  1.0: Want to plot the outline of a bw image;

function [] = plotBWOutline_v1_0(bw,style_text,lineWidth)

if(~exist('style_text','var'))
    style_text = '-';
end
if(~exist('lineWidth','var'))
    lineWidth = 0.25;
end
B = bwboundaries(bw);

hold on;
for i = 1:length(B)
    b = B{i};
    plot(b(:,2),b(:,1),style_text,'LineWidth',lineWidth);
end
hold off;