function ang = find_rot_angle_mod2(bodyTail_bw,min_con_reg) 
% input is rough bw mask of fish body (with or without tail fin - doesnt matter)
% output is ang - how much images needs to be rotated so that fish is approx aligned with horizontal 
%% Version History
%  mod1: using ang2 (nose to tail) angle only
%  mod2: using different method to calculate tail end

% erode body (find 'core' of body and tail - get rid of tail fin)
bodyTail_bw = bwmorph(bodyTail_bw,'erode',20); % not strictly nessesary
bodyTail_bw = bwareaopen(bodyTail_bw,0.25*min_con_reg);   % not strictly nessesary
 
% find nose and tail coor (for finding length of fish and later angle)
[y,x]=find(bodyTail_bw,10,'first');
nose_coor = [round(mean(x)),round(mean(y))];
% [y,x]=find(bodyTail_bw,1000,'last');
tail_coor = getBodyEndPoint_v1_0(bodyTail_bw,30,130);

% fish_length = tail_coor(1)-nose_coor(1);

% % point in the middle of the fish body at 0.35 of the full fish
% % length from the 'nose'
% midbody_coor(1) = nose_coor(1) + round(0.35*fish_length);
% 
% if(isnan(midbody_coor))
%     error('Orientation not detectable; likely need to adjust bwareaopen cutoff.');
% end
% First = find(bodyTail_bw(:,midbody_coor(1)),1,'first');
% Last  = find(bodyTail_bw(:,midbody_coor(1)),1,'last');
% 
% midbody_coor(2) = First+round(abs(First-Last)/2);
% 
% 
% % find angle of fish with respect to horizontal using midbody coor
% a = midbody_coor(1)-nose_coor(1);
% b = midbody_coor(2)-nose_coor(2);
% c = sqrt(a^2+b^2);
% ang1 = 0;
% if b>=0
%     ang1 = acos(a/c)*(360/(2*pi));
% elseif b<0
%     ang1 = -acos(a/c)*(360/(2*pi));
% end

% find angle of fish with respect to horizontal using tail coor
a = tail_coor(1)-nose_coor(1);
b = tail_coor(2)-nose_coor(2);
c = sqrt(a^2+b^2);
ang2 = 0;
if b>=0
    ang2 = acos(a/c)*(360/(2*pi));
elseif b<0
    ang2 = -acos(a/c)*(360/(2*pi));
end

% seems ang2 works best
ang = ang2;


% % plot midbody and tail coor
% figure(2)
% imshow(bodyTail_bw)
% hold on
% % plot([nose_coor(1) midbody_coor(1)],[nose_coor(2) midbody_coor(2)],'r-d')
% plot([nose_coor(1) tail_coor(1)],[nose_coor(2) tail_coor(2)],'g-o')
% pause(0.1);


end