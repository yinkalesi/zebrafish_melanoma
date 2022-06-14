%% findConfig_v1_3
%  Version 1.3
%  Author: Adeyinka Lesi
%  Date: 8/1/17

function [conf] = findConfig_v1_3(pris,counter,dist,max_depth,prev_nlocs)
% finds the unique configuration that minimizes the sum of the measure provided in
% dist. Useful for matching one set of objects to another
%% Version History
%  1.1: The old version sometimes got stuck as non-optimal configs were
%  searched. I am rewriting the code so it searches layer by layer - it
%  trys to make single changes to lower nflic instead of trying to minimize
%  nflic by going in-depth on the first conflict found...
%  1.2: Avoids exploring down 'bad' branches by making sure number of
%  conflict locations decreases - should run faster for complex systems
%  1.3: Looking for faster solution - can look at left-over solutions to
%  try and remove conflicts

[nr,nc] = size(pris);
depth = sum(counter-1); % number of changes made from original configuration

% MAX_DEPTH limits the recursive calls possible. The maximum number of
% possible calls should be O(nr^MAX_DEPTH)
if(~exist('max_depth','var'))
    max_depth = nc;
end

% if(nr>12)
%     MAX_DEPTH = min(4,nc);
%     fprintf('. ');
% %     disp(['Large number of recursion calls possible. nr = ' num2str(nr)]);
% %     disp(counter);
% end

% check that this configuration problem has a solutions
if(nr > nc)
    error('Unsolvable configuration problem');
end
% check to see if first configuration works
conf = pris(sub2ind(size(pris),1:nr,counter));
conflicts = cell(1,floor(0.5*nr));
nflic = 0;
nlocs = 0;
min_changes = 0;
% array to mark selected positions
repped_pos = zeros(1,nc);
flicted_pos = zeros(1,nc);
col_list = sort(pris(1,:));
for i = 1:length(col_list)
    check = find(conf==col_list(i));
    if(length(check)>=1)
        repped_pos(i) = 1;
        if(length(check)>1)
            % this is a conflict
            nflic = nflic+1;
            flicted_pos(i) = 1;
            nlocs = nlocs+length(check);
            min_changes = min_changes+length(check)-1;
            conflicts{nflic} = check;
        end
    end
end
conflicts = conflicts(1:nflic);
unrepped = find(repped_pos==0);
flicted = find(flicted_pos==1);
unassigned_i = sort([unrepped flicted]);
unassigned = col_list(unassigned_i);
nunas = length(unassigned);

ismain = 0; % variable differentiating first call of function from rest
if(~exist('prev_nlocs','var'))
    prev_nlocs = nlocs+1; % makes sure condition below is satisfied
    
    % if the condition above is satisfied, we are in the first function of
    % a series of recursive functions:
    ismain = 1;
end

% store all candidates for change in one array (locs)
locs = zeros(1,nlocs);
loci = 0;
for c = 1:nflic
    locs(1+loci:loci+length(conflicts{c})) = conflicts{c};
    loci = loci + length(conflicts{c});
end
locs = sort(locs);


if(ismain)
    depth_of_search = nc; 
    fprintf(['Starting %s. %i out of %i candidates are unrepresented. ' ...
        'Searching %i locations out of %i total.\n'],mfilename,...
        length(unrepped),nc,nlocs,nr);
else
    depth_of_search = max_depth;
    fprintf(['Starting %s. %i out of %i positions are unrepresented. ' ...
        'Searching %i locations.\n'],mfilename,length(unrepped),nc,...
        nlocs);
end

ncand = 0;
candidates = zeros(nr,nlocs);

if(depth<max_depth && prev_nlocs>nlocs)
    if(nflic>0)
        % check unrepped to see if unique locations can be found for each
        % position in a way that resolves conflicts
        if(nlocs < nr)
            % can use a subset matrix to examine only unassigned positions
            % can make a subpris matrix with locs as row and col representing
            % only unassigned positions
            subpris = zeros(nlocs,nunas);
            for ll = 1:nlocs
                [~,srow_i] = sort(pris(locs(ll),:));
                unas_i = srow_i(unassigned_i);
                [~,unas_ord] = sort(unas_i);
                subpris(ll,:) = unassigned(unas_ord);
            end
            subconf = findConfig_v1_3(subpris,ones(1,nlocs),...
                dist(locs,:),nunas,nlocs+1); % +1 allows more leeway
            if(~isempty(subconf))
                % change conf
                new_conf = conf;
                new_conf(locs) = subconf;
                % next candidate
                ncand = ncand+1;
                candidates(:,ncand) = new_conf;
            end
        else
            % no sub-matrix can be made
            % propose a change for each conflict location to find candidates
            for ll = 1:nlocs
                loc = locs(ll);
                new_counter = counter;
                new_counter(loc) = new_counter(loc)+1;
                % new_counter means an altered configuration is tested
                % the locs with number smaller than the present one are removed to
                % prevent redundant checking of configurations
                % depth_of_search prevents the algorithm from searching all the way
                % down every possible time - only one of the branches needs to
                % produce a candidate by depth_of_search
%                 disp(depth);
%                 disp(counter);
%                 disp(new_counter);
                new_conf = findConfig_v1_3(pris,new_counter,...
                    dist,depth_of_search,nlocs);
                if(~isempty(new_conf))
                    % store this viable configuration
                    ncand = ncand+1;
                    candidates(:,ncand) = new_conf;
                end
            end
        end
    else
        ncand = ncand+1;
        candidates(:,ncand) = conf;
    end
end

if(ncand==0)
    conf = [];
    if(ismain)
        error('Failed to find viable configuration');
    end
else
    % need selection criteria (should be minimum distance)
    totdist = zeros(1,ncand);
    for c = 1:ncand
        tconf = candidates(:,c);
        for j = 1:nr
            totdist(c) = totdist(c) + dist(j,tconf(j));
        end
    end
    [~,cc] = min(totdist);
    conf = candidates(:,cc);
end