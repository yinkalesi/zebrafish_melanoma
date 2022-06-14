function [bf] = MergeBF(bf1,bf2)

 % Merge bf1-2 take darkest pixels from each
    bf=bf1;
    [r,c]=size(bf);
    
        for kk=1:r
           for nn=1:c
               if bf1(kk,nn)<=bf2(kk,nn)
               bf(kk,nn) = bf1(kk,nn);
               elseif bf1(kk,nn)>bf2(kk,nn)
               bf(kk,nn) = bf2(kk,nn);
               end
           end
        end

end