function [x_coor,y_coor,x_coor_base,y_coor_base,x_coor_slope] = return_spine_coor_T_mod2(objs,num_points)
% mod2 - instead of having the same x-axis position, the reference points
% for different bodies have different positions based on the distance
% along the spinal cord - we want the transformed fish to have the same
% spinal cord positions (this should result in better transformations along
% the x-axis)

start = 0;
fin = Inf;
% mean_length = 0;
% mean_nose = 0;
mean_end = 0;
spine_lines_avg = zeros(length(objs),size(objs(1).spine,2));
arc_lengths = spine_lines_avg;
spine_slope_avg = spine_lines_avg;
for t = 1:length(objs)
    init_pos = 200;
    spine1 = sum(objs(t).spine(:,init_pos:end));
    spine1 = spine1(spine1~=0);
    
    max_spine = max(spine1);
    width1 = graythresh(double(spine1)/max_spine)*max_spine;
%     width1 = spine1(round(0.82*length(spine1)))*1.5;
    buffer = 20;
    start1 = find(spine1>width1,1,'last')+init_pos+buffer; % -1+1
    slut1 = find(spine1,1,'last')+init_pos-1-buffer;
    
    first_pos = objs(t).input_points(1,1);
    end_pos = objs(t).input_points(2,1);
    mean_end = mean_end + end_pos;
%     body_length = end_pos-nose_pos+1;
%     mean_nose = mean_nose + nose_pos;
%     mean_length = mean_length + body_length;
    start = max(start,start1);
    fin = min(fin,slut1);
    spine_line = zeros(1,size(spine_lines_avg,2));
    for i = start1-buffer:slut1+buffer
        spine_line(i) = getYPosition(objs(t).spine,i);
        if(spine_line(i)==0)
            % best guess is previous number
            spine_line(i) = spine_line(i-1);
        end
    end
    % smoothen curve
    avg_wid = 3; % this should definitely be smaller than buffer
    for i = start1-buffer:slut1+buffer
        low = max(i-avg_wid,start1-buffer);
        high = min(i+avg_wid,slut1+buffer);
        spine_lines_avg(t,i) = mean(spine_line(low:high));
    end
    spine_lines_avg(t,round(first_pos)) = objs(t).input_points(1,2);
    spine_lines_avg(t,round(end_pos)) = objs(t).input_points(2,2);
    
    npts = find(spine_lines_avg(t,:)~=0);
    spine_slope = zeros(1,size(spine_lines_avg,2));
    for ni = 1:length(npts)
        n = length(npts)-ni+1; % go backwards since zero is at tail
        % nn goes from length(npts) to 1
        if(n==length(npts))
            x = npts(n-2:n);
        elseif(n==1)
            x = npts(n:n+2);
        else
            x = npts(n-1:n+1);
        end
        del1 = x(2)-x(1);
        del2 = x(3)-x(2);
        del3 = x(3)-x(1);
        denom = del3*del2*del1;
        Ac1 = 2/denom;
        Ac2 = -(x(2)+x(3))/denom;
        Bc2 = -1/del2;
        Bc1 = (x(3)-x(1))*Bc2;
        if(n==1)
            range = x(1):x(2);
        else
            range = x(2):x(3);
        end
        for nni = range
            nn = range(end)-nni+range(1);
            A = nn*Ac1 + Ac2;
            B = A*Bc1 + Bc2;
            C = -A-B;
            dydx = [A B C]*spine_lines_avg(t,x)';
            spine_slope(nn) = dydx;
        end
    end
    % average out spine_slope to make a smooth curve
    avg_wid2 = 20;
    for r = npts(1):npts(end)
        low = max(r-avg_wid2,npts(1));
        high = min(r+avg_wid2,npts(end));
        spine_slope_avg(t,r) = mean(spine_slope(low:high));
    end
    % use averaged slope to get arc_length
    for ni = 1:npts(end)
        % go backwards
        n = npts(end)-ni+1;
        dydx = spine_slope_avg(t,n);
        arc_lengths(t,n) = arc_lengths(t,n+1) + sqrt(1+dydx^2);
    end
end


% mean_length = mean_length/length(objs);
% mean_nose = mean_nose/length(objs);
mean_end = mean_end/length(objs);
toflip = linspace(start,fin,num_points);
arc_coor = mean_end-toflip;

x_coor = zeros(length(objs),length(arc_coor));
x_coor_slope = x_coor;
for t = 1:length(objs)
    i_count = 0;
    first_pos = objs(t).input_points(1,1);
    end_pos = objs(t).input_points(2,1);
    for n = first_pos:end_pos
        if(i_count==length(arc_coor))
            break;
        end
        if(arc_lengths(t,n)<=arc_coor(i_count+1))
            i_count = i_count+1;
            x_coor(t,i_count) = n;
            x_coor_slope(t,i_count) = spine_slope_avg(t,n);
        end
    end
end     
y_coor = zeros(length(objs),length(arc_coor));
x_coor_base = zeros(1,length(arc_coor));
y_coor_base = zeros(1,length(arc_coor));
for t = 1:length(objs)
    for ni=1:length(arc_coor)
        y_coor(t,ni) = getYPosition(objs(t).spine,x_coor(t,ni));
    end
    x_coor_base = x_coor_base + x_coor(t,:);
    y_coor_base = y_coor_base + y_coor(t,:);
end
x_coor_base = round(x_coor_base/length(objs));         
y_coor_base = round(y_coor_base/length(objs));
         

end

function [ypos] = getYPosition(spine,xpos)
line = spine(:,xpos);
ind1 = find(line,1,'first');
ind2 = find(line,1,'last');
if(~isempty(ind1))
   ypos = 0.5*(ind1+ind2);
else
   ypos = 0;
end
end