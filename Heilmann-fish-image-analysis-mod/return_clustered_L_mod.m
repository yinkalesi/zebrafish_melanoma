function L_clustered = return_clustered_L_mod(BW,distMatrix,threshold,time)
% returnClusteredL() takes black and white image with corresponding
% distmatrix and determines which objects are closer than threshold pixels
% and labels these blops with the same index in L_clustered

if isempty(distMatrix)==0 && isempty(BW)==0
    % distMatrix and BW was not empty so bw need to be clustered

    L_clustered = bwlabel(BW);
    distVec     = squareform(distMatrix); % returns a vector form of distmatrix (which linkage needs)
    LIN         = linkage(distVec); % returns a matrix LIN that encodes a tree of hierarchical clusters of the rows of distMatrix (row/col number of object is the same as index of object in Labeled matrix)
    % LIN has size (numBlobs-1)x3. 1st and 2nd column hold index of two
    % dofferent objects and 3rd column has the pixels distance/threshold
    % which would make these two objects clustered

    numBlobs = max(L_clustered(:));
    
    % Iterate through 3rd coulumn of LIN and compare with threshold to find
    % out which objects to cluster
    countJoins=0; % count number of 'joins' we need to make (max possible will be numBlops-1)
    for g = 1:length(LIN(:,3))
        if LIN(g,3)<threshold
            
            countJoins = countJoins+1;
            
            L_clustered(L_clustered==LIN(g,2)) = numBlobs+countJoins; % replace different index in these two objects with the same number
            L_clustered(L_clustered==LIN(g,1)) = numBlobs+countJoins; % note - this number starts at numBlops + 1 so it will be different from other index already there
        end
    end
    
    % We want the largest index in L_clustered to be equal to the number of
    % clustered blobs (number of objects=tumors/metastases)
    % so we need to go back and replace the large index numbers (we got above) with the
    % lower index numbers which are now freed up!
    for i=1:numBlobs-countJoins % the new number of labeled clutered regions is: numBlobs-countJoins
        if isempty(find(L_clustered==i,1))==1 % if the number i is not being used as label in L_clustered then take the largest index used and replace it with i
            L_clustered(L_clustered==max(L_clustered(:))) = i;
        end    
    end
    
elseif isempty(distMatrix)==1 && isempty(BW)==0 % if distMatrix is empty just return normal labeled matrix

    L_clustered = bwlabel(BW);
    
elseif isempty(BW)==1 % if BW is an empty matrix - return empty matrix
    
    L_clustered = [];
    
end

% if time==1 
%     L_clustered = BW;% cluster everything if time is day one (we do not want to count a fragmented implant as several objects...)
% end

end % end of function

