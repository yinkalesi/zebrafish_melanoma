
function fish_length = return_fish_length(body)


        % find nose coor so that images can be cropped
        [y,x]=find(body,10,'first');
        nose_coor = [mean(x),mean(y)];
        [y,x]=find(body,10,'last');
        tail_coor = [mean(x),mean(y)];
        
        fish_length = tail_coor(1)-nose_coor(1);


end