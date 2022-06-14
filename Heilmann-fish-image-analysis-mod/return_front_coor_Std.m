function [input_points,base_points] = return_front_coor_Std(body,body_std)


body_overlap = body;
body_overlap(body_std==0)=0;

fish_length = return_fish_length(body_overlap);

[~,start]=find(body_overlap,1,'first');

slut = round(fish_length*0.98);

x_coor = round(linspace(start+20,slut-10,20));


yup_coor = zeros(1,length(x_coor));

ymiddle_coor = zeros(1,length(x_coor));

ydown_coor = zeros(1,length(x_coor));

yupStd_coor = zeros(1,length(x_coor));

ymiddleStd_coor = zeros(1,length(x_coor));

ydownStd_coor = zeros(1,length(x_coor));


for i=1:length(x_coor)
        yup_coor(i) = find(body(:,x_coor(i)),1,'first');
        yupStd_coor(i) = find(body_std(:,x_coor(i)),1,'first');

        ymiddle_coor(i) = 400;
        ymiddleStd_coor(i) = 400;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO SOMETHING

        
        ydown_coor(i) = find(body(:,x_coor(i)),1,'last');
        ydownStd_coor(i) = find(body_std(:,x_coor(i)),1,'last');    
end

x_coor = [x_coor x_coor x_coor];

y_coor = [yup_coor ymiddle_coor ydown_coor];


ydownStd_coor(12:17) = max([ydownStd_coor(12:17);ydown_coor(12:17)]); % use the fish own body shape (not the standard shape) in the tumor area to avoid huge deformations of tumor region 

%yupStd_coor(14:16) = max([yupStd_coor(14:16);yup_coor(14:16)]); % use the fish own body shape (not the standard shape) near the back fin area to avoid huge deformations


y_coor_base = [yupStd_coor ymiddleStd_coor ydownStd_coor];


input_points = zeros(length(x_coor),2);
base_points   = zeros(length(x_coor),2);

input_points(:,1) = x_coor;
input_points(:,2) = y_coor;

base_points(:,1) = x_coor;
base_points(:,2) = y_coor_base;


 
end

