function y_coor = return_Y_coor(x_coor,spine_bw)



         x = spine_bw(:,x_coor);
         s = circshift(x,1);
         y = abs(x-s); % ones at boundary
         
         if sum(y)==2 % only two boundaries - first and last is good enough
             y_coor=round((find(x,1,'last')-find(x,1,'first'))/2+find(x,1,'first')); 
             
             hej=2;
         elseif sum(y)==4
             hej=4;

             [z1,z2]=find(y); % z2 now hold indices of all boundary points
             z3 = circshift(z2',1)';
             c=z2-z3;
             c=c(2:2:4);
             [m,i]=max(c); % i is 1 if first stretch of ones is biggest and 2 if second strech of ones is biggest
             temp = z1(i*2-1:i*2);% get the right boundary points from z2
             
             y_coor = mean(temp); % mean should give point in between
         else
             sum(y)
            display('too many boundary points cannot find spine control point!!!!');
             
            y_coor =0;
            
         end
         
         y_coor;
end