function [x_coor,y_coor1,y_coor_base] = return_spine_coor_Std(spine_bw,spine_bw_std)

     spine   = sum(spine_bw(:,90:end));
     spineStd = sum(spine_bw_std(:,90:end));

     % get rid of zeros
     index=find(spine);
     spine = spine(index);
     indexStd=find(spineStd);
     spineStd = spineStd(indexStd);
              
         
     width   = spine(round(0.9*length(spine)))*1.5;
     widthStd = spineStd(round(0.9*length(spineStd)))*1.5;
     
     
     start   = find(spine<width,1) +90;
     startStd = find(spineStd<widthStd,1) +90;

     slut = find(spine,1,'last')+90;
     slutStd = find(spineStd,1,'last')+90;
     
     start = max([start startStd]) ;
     slut = min([slut slutStd]);
     
     x_coor = round(linspace(start+1,slut-1,5));
     
     y_coor1     = x_coor.*0;
     y_coor_base = y_coor1;


         for n=1:length(x_coor)
             y_coor1(n)     = return_Y_coor(x_coor(n),spine_bw);             
             y_coor_base(n) = return_Y_coor(x_coor(n),spine_bw_std);

         end
         

end