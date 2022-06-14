function [gfp] = Merge(gfp1,gfp2)

 % Merge gfp1-2 take briftest pixels from each
    gfp=gfp1;
    [r,c]=size(gfp);
    
        for kk=1:r
           for nn=1:c
               if gfp1(kk,nn)>=gfp2(kk,nn)
               gfp(kk,nn) = gfp1(kk,nn);
               elseif gfp1(kk,nn)<gfp2(kk,nn)
               gfp(kk,nn) = gfp2(kk,nn);
               end
           end
        end

end