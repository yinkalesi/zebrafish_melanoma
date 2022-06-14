function V_labeled = return_labeled_voronoi_mod2(x,y,rows,cols,threshold)
% returnLabeledVoronoi takes x- and y-coordinates of N>2 points and makes a
% matrix of size (rows,cols) with the different voronoi regions (defined by
% the N points) labelled 1,2,3 ... N.
%
% each region is seperated from the other by a line of 0's of width 1-2 pixels (which is why the input points should not be too close!)
% Detailed explanation goes here


% TEST IF EVERYTHING IS IN ORDER
% test whether some points are too close to do voronoi!
for ii=1:length(x)
    for jj=1:length(x)
        dist = sqrt((x(ii)-x(jj))^2 +(y(ii)-y(jj))^2);
        if ii~=jj && dist <threshold
            % This happens because different objects can be arbitrarily
            % close. It doesn't necessarily cause a big problem so I will
            % remove this message (any error should be caught later)
%             disp('some points are too close to use returnLabeledVoronoi!')
        end
    end
end

% test if there is enough points to do voronoi!
if length(x)<2
    % this shouldn't happen
    error('Incorrect point specification for voronoi');
elseif length(x)==2
    warning('Not enough points to use voronoi; using alternative method');
    % this means there are two points; need to find the line representing
    % the plane between the points
   
    % need a point on separation line - this will be set to be the point on
    % the connecting line half way between voronoi points
    xp = mean(x);
    yp = mean(y);
    % verticies of line should be where the line intersects with edge of
    % image
    if((x(2)-x(1)) == 0)
        vx = [1 cols]';
        vy = [yp yp];
    elseif((y(2)-y(1)) == 0)
        vy = [1 rows]';
        vx = [xp xp];
    else
        % slope of separation line is perpendicular to connecting line
        slope = -1*(x(2)-x(1))/(y(2)-y(1));
        vy = [1 rows]';
        vx = xp + (vy-yp)/slope;
    end
else
    
    % MATLABs voronoi function returns the finite vertices of the Voronoi edges
    % in vx and vy (coordinates of end points of lines between regions)
    [vx,vy]=voronoi(x,y);
end
% Initialize empty matrix
V = zeros(rows,cols);

% We want the outer edges of the voronoi lines to reach the edge of the
% image. This means we should identify which nodes are boundary points
% and extend the line segments belonging to those nodes
vbounds = boundary([vx(:); x],[vy(:); y]);

% will check each line segment to see which has a node on the convex hull
% and those line segments will be extended

% iterate through all voronoi edges 
for h=1:size(vx,2)
    % check if nodes are corner of convex all
    iscorner1 = any(ismember(vbounds,2*h-1));
    iscorner2 = any(ismember(vbounds,2*h));
    
    teller = vy(1,h)-vy(2,h);
    nevner = vx(1,h)-vx(2,h);
    
    if nevner==0
        % if slope is Inf (nevner ==0, edge is vertical) we need to generate xx and yy in a
        % diferent way
        
        % change YMIN and YMAX based on corner status
        % if the node of a line segment is a corner, then that node needs to be
        % pushed to the edge of the image (YMIN=1 or YMAX=rows respectively)
        if(vy(1,h) < vy(2,h))
            if(iscorner1)
                YMIN = 1;
            else
                YMIN = vy(1,h);
            end
            if(iscorner2)
                YMAX = rows;
            else
                YMAX = vy(2,h);
            end
        else
            if(iscorner1)
                YMAX = rows;
            else
                YMAX = vy(1,h);
            end
            if(iscorner2)
                YMIN = 1;
            else
                YMIN = vy(2,h);
            end
        end
        
        yy = YMIN:0.5:YMAX;
        xx = vx(1,h)*ones(1,length(yy));
    else        
        % slope of voronoi edge
        slope = teller/nevner;
        
        % change XMIN and XMAX based on corner status
        % if the node of a line segment is a corner, then that node needs to be
        % pushed to the edge of the image (XMIN=1 or XMAX=cols respectively)
        if(vx(1,h) < vx(2,h))
            if(iscorner1)
                XMIN = 1;
            else
                XMIN = vx(1,h);
            end
            if(iscorner2)
                XMAX = cols;
            else
                XMAX = vx(2,h);
            end
        else
            if(iscorner1)
                XMAX = cols;
            else
                XMAX = vx(1,h);
            end
            if(iscorner2)
                XMIN = 1;
            else
                XMIN = vx(2,h);
            end
        end
    
        % We need to draw a line such that the spacing is less than 1 (since
        % we want to paint each pixel along the way). Therefore, dx =
        % min(0.5,0.5/|slope|) and dy = min(0.5,0.5*|slope|)
        dx = min(0.5,0.5/abs(slope));
        
        % x-coor span for drawing edge/line on V
        xx = XMIN:dx:XMAX;
        
        % intersection with y-axis of edge/line
        b = vy(1,h)-vx(1,h).*slope;
        
        % y - values of line
        yy = xx.*slope+b;
    end

    % get rid of xx's and yy's that are outside V. (Set them to 1) ... some
    % of these might be overkill... hmmm
    xx(yy>rows)=1;
    yy(yy>rows)=1;
    xx(yy<1)=1;
    yy(yy<1)=1;
    xx(xx>cols)=1;
    yy(xx>cols)=1;
    xx(xx<1)=1;
    yy(xx<1)=1;
    
    % make new temp matrix with 1's (same size as V)
    temp = zeros(rows,cols);

    % draw edge/line with 1's in temp matrix of 0's   
    for i=1:length(xx)
        temp(round(yy(i)),round(xx(i)))=1;
    end
    % remove the white pixel set to 1 by the [xx,yy]=[0,0]
    temp(1,1)=0; 

    % make line fatter by dilating
    temp=bwmorph(temp,'dilate',1);
    
    % add line to V
    V = V + temp;
    
end

% set everything > 0 to 1
V(V>0)=1;
% make negative
V=V.*-1+1;

% make lines thinner by thickening regions
V = bwmorph(V,'thicken',5);
% label each region with 1,2,3,....,N      
V_labeled=bwlabel(V);

end

