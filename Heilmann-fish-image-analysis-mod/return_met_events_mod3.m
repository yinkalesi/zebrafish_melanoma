function [cmmets,index] = return_met_events_mod3(objs)
%returnMetEvents returns list with center of mass of all metastatic events which took place between two time points 
% cmD1 list of center of mass of mets/tumor which arrose between day 0 and day1 (including the injection site!!)
% cmD7 list of center of mass of mets which arrose between day 1 and day 7 
% cmD14 list of center of mass of mets which arrose between day 7 and day 14 

%% Version History
%  mod3 is same as mod2 for now

nt = length(objs);

nmets = zeros(1,nt);

index = cell(1,nt);

cmmets = cell(1,nt);

for rt=1:nt
    t = nt-rt+1;
    indexList = nonzeros(unique(objs(t).L_t1));
    % go through objects on day 14
    for i=1:length(indexList)
        
        % make image with only object i ...
        L_i = return_sub_listL(objs(t).L_t1,indexList(i));
        
        % look back in time to see if any previous time has an object
        % that matches
        has_match = 0;
        
        if(t>1)
            % ...(then dilate it so change of overlap is increased)
            L_i = bwmorph(L_i,'dilate',0.5*objs(t).threshold);
            for tt = 1:(t-1)
                tempL = objs(t-tt).L_t1;
                tempL(L_i==0)=0;
                % get list of indexes for objects at day 7 contributing to object i at day 14
                if(~isempty(nonzeros(unique(tempL))))
                    has_match = 1;
                end
            end
        end
           
        if ~has_match % if there are no objects within the dilated object #i of D14, at day 7 it means we had a metastatic event between day 7 and day 14
            if isempty(find(index{t}==indexList(i),1))==1 % only add metastatic event to list if it is not there already!
                nmets(t) = nmets(t) + 1; % count metastatic event
                L_i = double(L_i);
                L_i(L_i>0)=2;
                r=regionprops(L_i,'Centroid');
                cmmets{t}(nmets(t),:) = r(2).Centroid; % r(2) because you just set number>0 to 2's!!
                index{t}(nmets(t))=indexList(i);
            end
        end
    end
end

end

