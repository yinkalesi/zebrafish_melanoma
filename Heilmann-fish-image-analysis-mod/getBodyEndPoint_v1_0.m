function [endpoint,orientation] = getBodyEndPoint_v1_0(body,shift1,shift2)
%% getTailEndPoint
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 9/1/17
%  Project: Fish Image Analysis
%  find the endpoint of the body (the center of the end plane)

[lastys,lastxs] = find(body,100,'last');
if(isempty(lastxs))
    error('Body image is blank');
end
lastx = mean(lastxs);
lasty = mean(lastys);

% approximate tail end of body with two lines (the end of body is curved so
% will try and avoid this region);
npts = 10;
% shift1 = 50;
% shift2 = 150;
xpos = round(linspace(lastx-shift2,lastx-shift1,npts));
yups = zeros(1,npts);
ydowns = zeros(1,npts);

for i = 1:npts
    line = body(:,xpos(i));
    yups(i) = find(line,1,'first');
    ydowns(i) = find(line,1,'last');
end

pup = polyfit(xpos,yups,1);
pdown = polyfit(xpos,ydowns,1);

% get orientation
avg_slope = 0.5*(pup(1)+pdown(1));
perp_slope = -1/avg_slope;
orientation = atan(perp_slope); % perpendicular angle

% you can create a box using yup, ydown and the perpendicular line that
% goes though lastx. The endpoint is the center of the right edge of the box

% find intersects (box corners)
% assume yup and ydown are not vertical
if(isinf(pup(1)) || isinf(pdown(1)))
    error('Body lines are vertical');
end
% if perpendicular line is vertical (avg_slope = 0), need alternative
if(avg_slope == 0)
    isc1x = lastx;
    isc2x = lastx;
else
    isc1x = (lasty - pup(2) - lastx*perp_slope)/(pup(1)-perp_slope);
    isc2x = (lasty - pdown(2) - lastx*perp_slope)/(pdown(1)-perp_slope);
end
isc1y = polyval(pup,isc1x);
isc2y = polyval(pdown,isc2x);

% endPoint is halfway up the line
halfdist = 0.5*sqrt((isc1x-isc2x)^2+(isc1y-isc2y)^2);
endpoint = round([isc1x+sign(isc2x-isc1x)*cos(orientation)*halfdist, ...
    isc1y+abs(sin(orientation))*halfdist]);

% figure(4);
% imshow(body);
% hold on;
% plot(xpos,yups,'bo');
% plot(xpos,ydowns,'ro');
% plot([xpos isc1x],polyval(pup,[xpos isc1x]),'b-');
% plot([xpos isc2x],polyval(pdown,[xpos isc2x]),'r-');
% plot([isc1x endpoint(1) isc2x],[isc1y endpoint(2) isc2y],'yo-');
% plot(lastx,lasty,'*');
end

