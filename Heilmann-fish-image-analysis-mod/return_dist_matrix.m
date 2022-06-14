function distMatrix = return_dist_matrix(BW)
% returnDistMatrix(BW) takes a black/white (0/1) matrix/image and return a pairwise distance
% matrix holding the shortest distances between all seperate regions with
% 1's. (distMatrix is a symmetric matrix)
%
% Distances are shortest possible distances between points on the
% boundary of both objects
% if BW is all 0's returnDistMatrix(BW) will return an empty matrix []


% get boundary points of the blops in BW
B = bwboundaries(BW,'noholes');

% Only contruct a distance matrix if there is more than one blop (aka seperate region with 1's)
if length(B)>1 % length(B) is number of blops
    
    % initialize
    distMatrix = zeros(length(B));
    
    for m = 1:length(B) % iterate through all different blops
        for k = m:length(B) % iterate through blops that have not allready apeared in the above iteration
            
            if m~=k % only find min dist between different objects
                
                object1=B{m}; % holds boundary points of object #1
                
                object2=B{k}; % holds boundary points of object #2
                
                % put all x-coor in first row - all y-coor in 2 row and put
                % 1 and 2 in 3rd row to keep track of which point came from
                % what object
                points = [object1(:,1)' object2(:,1)' ; object1(:,2)' object2(:,2)' ; ones(1,length(object1)) 2*ones(1,length(object2))];
                points = points';
                
                Pd = pdist(points(:,1:2)); % returns pairwise distance between all points in vector format
                Pd_sq = squareform(Pd);    % converts vector to matrix format Pd_sq(i,j) gives dist between point i and j
                
                [row,col] = size(Pd_sq);
                
                % Start out with rediculously large number - replace
                % everytime you find a smalle dist (looking for the smallest possible)
                mindist = 1e10;
                
                for i=1:row
                    for j=1:col
                        if mindist>Pd_sq(i,j) && points(i,3)~=points(j,3) % only considder distances between points from different objects!!!
                            mindist=Pd_sq(i,j);
                        end
                    end
                end
                
                distMatrix(k,m) = mindist; % save smallest distance in distMatrix
            end
        end
    end
    
    % fill both upper and lower triangle
    temp = distMatrix';
    distMatrix = distMatrix + temp;
else
    
    % Make sure the function have something called 'distMatrix' to return
    distMatrix = [];

end


end


