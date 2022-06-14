function [x0,place_dist,place_rad,dist0,prim_list,meta_list,inferred_data,size_projections] = getOriginSizes_v1_0(key,data,prim_acceptance_thresh)
%% getOriginSizes_v1_0
%  version 1.0
%  Author: Adeyinka Lesi
%  Date: 5/15/20
% get TIME_ZERO_SIZES using RK4 (step back in time using given velocities)
% need to handle tumors that are known to be present but are below
% detection limit. That implies that the specified size in x_init is
% irrelevant for these 'invisible' tumors and they will be placed based on
% a putative distribution
%% Version History
%  1.0: from getTimeZeroSizes_v4_1, adding ability to track metastasis and
%  their times of origin
%  11/7/20: noticed issue when immunity starts early enough that it effects
%  trajectory calculation. getSizeTrajectory_AlwaysSmaller_v4_0 is set to
%  stop calculating at any point v(x)<0 so trajectory gives back the
%  initial size for each time point. Additionally there was an if statement
%  traj(end) < x_ti(j) for setting x0i. Since traj(end)==x_ti(j), x0i was
%  being set to 1. Changed to getSizeTrajectory_AlwaysSmaller_v4_1
%  2/3/21: adding extra comments
%  3/31/21: made it so place_dist/rad at time 1 is set based on values in
%  key.PLACEMENT_DIST/RADIUS


if(isfield(key,'QUIET_MODE') && ~key.QUIET_MODE)
    display(['> ' mfilename]);
end

% initialize: 
% prim_acceptance_thresh refers to how many time points the code
% will look at to check for inferrable primaries
if(~exist('prim_acceptance_thresh','var'))
    prim_acceptance_thresh = 3;
end

nt = length(data.t);

% this is a standard adjustment to the data distribution for cases when the
% primary isn't included at every time point
if(isfield(data,'cdf_adj'))
    p_adj = data.cdf_adj(1);
else
    p_adj = 0;
end
x_t = cell(1,nt); % tumor sizes at each
dist = cell(1,nt); % normalized tumor count at each time
cdf = cell(1,nt); % cumulative distribution
x_t{1} = data.x(data.selector{1});
dist{1} = data.dist(data.selector{1},1);
cdf{1} = p_adj+data.cdf(end,1)-[0; data.cdf(data.selector{1}(1:end-1),1)];
x0_size = length(x_t{1});
for n = 2:nt
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(n);
    else
        p_adj = 0;
    end
    x_t{n} = data.x(data.selector{n});
    dist{n} = data.dist(data.selector{n},n);
    cdf{n} = p_adj+data.cdf(end,n)-[0; data.cdf(data.selector{n}(1:end-1),n)];
    x0_size = max(x0_size,length(x_t{n}));
end

inferred_data = data;

temp_data = struct('x',cell(1,nt),'dist',cell(1,nt),'ccdf',...
    cell(1,nt),'step_ccdf',cell(1,nt),'labels',cell(1,nt),...
    'detect_time',cell(1,nt));

x0 = zeros(x0_size,1);
dist0 = zeros(x0_size,1);
prim_list = zeros(x0_size,3);
known_list = zeros(x0_size,5);
meta_list = zeros(x0_size,5);
backtrack_sizes = cell(1,nt);
origin_times = cell(1,nt);
nx0 = 0;
nkwn = 0;
nmet = 0;
place_dist = zeros(size(x0));
place_rad = zeros(size(x0));
i_selected = cell(1,nt);
x_t1_low = 0;
ntum_inc = 0; % a count of normalized number of tumors known to be included in the distribution
low_check1 = 0; % low limit for first time point
low_check2 = 1; % low limit for other time points
for i = 1:nt
    x_ti = x_t{i};
    x0i = zeros(length(x_ti),1);
    if(i==1)
        place_disti = key.PLACEMENT_DIST;
        place_radi = key.PLACEMENT_RADIUS;
    else
        place_disti = zeros(size(x0));
        place_radi = zeros(size(x0));
    end
    ts = data.t(i);
    prev_sizes_store = zeros(length(x_ti),i);
    birth_times = zeros(length(x_ti),1);
    if(ts>0)
        dt = ts/round(ts/key.TIME_STEP)/2;
        for j = 1:length(x_ti)
            if(x_ti(j)==0)
                % indication disappeared tumor
                x0i(j) = 0;
                prev_sizes = zeros(1,i);
                prev_sizes_store(j,:) = prev_sizes;
                birth_times(j) = -1;
            else
                [traj,ttraj] = getSizeTrajectory_AlwaysSmaller_v4_1(key,x_ti(j),-dt,ts);
                if(traj(end) <= x_ti(j))
                    x0i(j) = max(1,traj(end));
                else
                    % net negative velocity
                    % 1) parameters don't make sense
                    % 2) will set initial point to 1 (for now) - effect of this
                    % depends on low_checks
                    x0i(j) = 1;
                end
                % store predicted sizes at previous times
                prev_ti = length(traj)-round(data.t(1:i)/dt);
                prev_sizes = traj(prev_ti);
                prev_sizes(prev_sizes<1) = 0;
                prev_sizes_store(j,:) = prev_sizes;
                % get birth time
                birth_ind = find(traj>1,1,'last');
                if(~isempty(birth_ind))
                    birth_times(j) = ttraj(birth_ind);
                end
            end
        end
    else
        for j = 1:length(x_ti)
            x0i(j) = x_ti(j);
            prev_sizes = x_ti(j);
            prev_sizes_store(j,:) = prev_sizes;
            %             if(x_ti(j) >= lowx)
            %                 x0i(j) = x_ti(j);
            %             else
            %                 x0i(j) = 1+0.5*(lowx-1);
            %                 place_disti(j) = unk_dist_key.(unk_dist);
            %                 place_radi(j) = 0.5*(lowx-1);
            %             end
        end
    end
    
    % test to see what initial tumors should be accepted
    
    % get rid of x0 values < low_check1
    i_above = find(x0i>low_check1);
    if(isempty(i_above))
        fprintf('Using 0 initial tumor values from time %i\n',i);
    else
        if(i == 1)
            if(length(i_above) < length(x0i))
                fprintf('Using %i out of %i initial tumor values from time %i\n',length(i_above),length(x0i),i);
            else
                fprintf('Using all %i initial tumor values from time %i\n',length(x0i),i);
            end
            
            nai = length(i_above);
            x0(nx0+1:nx0+nai) = x0i(i_above);
            dist0(nx0+1:nx0+nai) = dist{i}(i_above);
            prim_list(nx0+1:nx0+nai,:) = [ones(nai,1)*i i_above x_ti(i_above)];
            place_dist(nx0+1:nx0+nai) = place_disti(i_above);
            place_rad(nx0+1:nx0+nai) = place_radi(i_above);
            i_selected{i} = i_above;
            nx0 = nx0+nai;
            
            % this is not needed...
            %             known_sizes = x0i(i_above);
            %             known_origin_time = birth_times(i_above);
            %             known_list(nkwn+1:nkwn+nai,:) = [ones(nai,1)*i i_above known_sizes known_origin_time];
            %             nkwn = nkwn+i_above;
            
            if(nx0 > 0)
                x_t1_low = x0(nx0); % lowest size at time 1
                ntum_inc = sum(dist{1}(i_above)); % normalized number of tumors included
            elseif(~isempty(x_t{1}))
                x_t1_low = x_t{1}(end);
                ntum_inc = 0;
            else
                x_t1_low = inf;
                ntum_inc = 0;
            end
            %             x_t1_low = inf;
            %             n_t1 = 0;
        else
            % want to select init values larger than 1 for the case when
            % there are excess tumors at later times
            
            if(i <= prim_acceptance_thresh)
                prim_sel = find((x0i>low_check2).*(cdf{i}>ntum_inc)); % labeled as new (formerly undetected) primary
                meta_sel = find((x0i>0).*(x0i<=low_check2).*(cdf{i}>ntum_inc)); % this data is attributed to new (formerly undetected) metastases
                known_sel = find((x0i>0).*(cdf{i}<=ntum_inc)); % this data is attributed to previously identified tumors (prim or met)
            else
                prim_sel = [];
                meta_sel = find((x0i>0).*(cdf{i}>ntum_inc)); % this data is attributed to new (formerly undetected) metastases
                known_sel = find((x0i>0).*(cdf{i}<=ntum_inc)); % this data is attributed to previously identified tumors (prim or met)
            end
            if(i <= prim_acceptance_thresh)
                if(length(prim_sel) < length(x0i))
                    fprintf('Using %i out of %i initial tumor values from time %i\n',length(prim_sel),length(x0i),i);
                else
                    fprintf('Using all %i initial tumor values from time %i\n',length(x0i),i);
                end
            end
            
            % update primaries
            nai = length(prim_sel);
            x0(nx0+1:nx0+nai) = x0i(prim_sel);
            dist0(nx0+1:nx0+nai) = dist{i}(prim_sel);
            prim_sizes = x_ti(prim_sel);
            prim_list(nx0+1:nx0+nai,:) = [ones(nai,1)*i prim_sel prim_sizes];
            place_dist(nx0+1:nx0+nai) = place_disti(prim_sel);
            place_rad(nx0+1:nx0+nai) = place_radi(prim_sel);
            i_selected{i} = prim_sel;
            nx0 = nx0+nai;
            
            % update known primaries and metastasis
            nki = length(known_sel);
            known_sizes = x_ti(known_sel);
            known_dist = dist{i}(known_sel);
            known_origin_time = birth_times(known_sel);
            known_list(nkwn+1:nkwn+nki,:) = [ones(nki,1)*i known_sel known_sizes known_dist known_origin_time];
            nkwn = nkwn+nki;
            
            % update metastases
            nmi = length(meta_sel);
            meta_sizes = x_ti(meta_sel);
            meta_dist = dist{i}(meta_sel);
            meta_origin_time = birth_times(meta_sel);
            meta_list(nmet+1:nmet+nmi,:) = [ones(nmi,1)*i meta_sel meta_sizes meta_dist meta_origin_time];
            nmet = nmet+nmi;
            
            % update total number of tumors in distribution
            ntum_inc = ntum_inc+sum(dist{i}(prim_sel))+sum(dist{i}(meta_sel));
        end
    end
    backtrack_sizes{i} = prev_sizes_store;
    origin_times{i} = birth_times;
end

if(nx0 == 0)
    x0(1) = 1;
    nx0 = nx0+1;
    warning('None of the initial tumor values were used. Using one size 1 tumor');
end

x0 = x0(1:nx0);
[~,ord] = sort(x0(1:nx0));
x0 = x0(ord);
dist0 = dist0(ord);
prim_list = prim_list(ord,:);
place_dist = place_dist(ord);
place_rad = place_rad(ord);
known_list = known_list(1:nkwn,:);
meta_list = meta_list(1:nmet,:);

size_projections = struct('backtrack_sizes',backtrack_sizes,'origin_times',origin_times);

% construct new data distribution (with inferred primaries and with and 
% without inferred meta
% add backtrack_sizes and origin_times to structure

% construct distributions
% 1) known (visible prim and meta)+inferred prim
% 2) known (visible prim and meta)+inferred prim+inferred meta
nxt = 0;
for i = 1:nt
    known_i = known_list(:,1)==i;
    prim_i = prim_list(:,1)==i;
    meta_i = meta_list(:,1)==i;
    infp_i = prim_list(:,1)>i;
    infm_i = (meta_list(:,1)>i)&(meta_list(:,end)<=data.t(i));
    n_known = sum(known_i);
    n_newprim = sum(prim_i);
    n_newmeta = sum(meta_i);
    n_infprim = sum(infp_i);
    n_infmeta = sum(infm_i);
    
    xi_lab = [ones(n_infmeta,1)*2; ones(n_infprim,1)*1; zeros(n_newmeta+n_newprim+n_known,1)];
    xi_orig = zeros(size(xi_lab));
    xi = zeros(size(xi_lab));
    disti = zeros(size(xi_lab));
    nxi = 0;
    xi(end-nxi-n_known+1:end-nxi) = backtrack_sizes{i}(known_list(known_i,2),i);
    disti(end-nxi-n_known+1:end-nxi) = dist{i}(known_list(known_i,2),1);
    nxi = nxi+n_known;
    xi(end-nxi-n_newprim+1:end-nxi) = backtrack_sizes{i}(prim_list(prim_i,2),i);
    disti(end-nxi-n_newprim+1:end-nxi) = dist{i}(prim_list(prim_i,2),1);
    nxi = nxi+n_newprim;
    xi(end-nxi-n_newmeta+1:end-nxi) = backtrack_sizes{i}(meta_list(meta_i,2),i);
    disti(end-nxi-n_newmeta+1:end-nxi) = dist{i}(meta_list(meta_i,2),1);
    nxi = nxi+n_newmeta;
    xi_orig(end-nxi+1:end) = ones(nxi,1)*i;
    
    infp_inds = find(infp_i);
    for j = 1:n_infprim
        xi(end-nxi) = backtrack_sizes{prim_list(infp_inds(j),1)}(prim_list(infp_inds(j),2),i);
        disti(end-nxi) = dist{prim_list(infp_inds(j),1)}(prim_list(infp_inds(j),2),1);
        xi_orig(end-nxi) = prim_list(infp_inds(j),1);
        nxi = nxi+1;
    end
    infm_inds = find(infm_i);
    for j = 1:n_infmeta
        xi(end-nxi) = backtrack_sizes{meta_list(infm_inds(j),1)}(meta_list(infm_inds(j),2),i);
        disti(end-nxi) = dist{meta_list(infm_inds(j),1)}(meta_list(infm_inds(j),2),1);
        xi_orig(end-nxi) = meta_list(infm_inds(j),1);
        nxi = nxi+1;
    end
    
    % sort by size
    [~,ord] = sort(xi);
    xi = xi(ord);
    disti = disti(ord);
    xi_lab = xi_lab(ord);
    xi_orig = xi_orig(ord);
    
    % calculate cdf...
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(i);
    else
        p_adj = 0;
    end
    cdfi = getStepIntegral(xi,disti);
    ccdfi = p_adj+cdfi(end)-[0; cdfi];
    ccdf2i = zeros(2*length(xi),3); % double cdf for making plots
    ccdf2i(2*(1:length(xi))-1,1) = xi;
    ccdf2i(2*(1:length(xi)),1) = xi;
    ccdf2i(2*(1:length(xi))-1,2) = ccdfi(1:end-1);
    ccdf2i(2*(1:length(xi)),2) = ccdfi(2:end);
    ccdf2i(2*(1:length(xi))-1,3) = xi_lab;
    ccdf2i(2*(1:length(xi)),3) = xi_lab;
    
    temp_data(i).x = xi;
    temp_data(i).dist = disti;
    temp_data(i).ccdf = ccdfi(1:end-1);
    temp_data(i).step_ccdf = ccdf2i;
    temp_data(i).labels = xi_lab;
    temp_data(i).detect_time = xi_orig;
    
    nxt = nxt+nxi;
    
    %     % plott
    %     figure;
    %     semilogx(x_t{i},cdf{i},'--ob',ccdf2i(:,1),ccdf2i(:,2),'-g',...
    %         ccdf2i(ccdf2i(:,3)==1,1),ccdf2i(ccdf2i(:,3)==1,2),'sr',...
    %         ccdf2i(ccdf2i(:,3)==2,1),ccdf2i(ccdf2i(:,3)==2,2),'dr');
    %     if(n_infprim>0 && n_infmeta>0)
    %         legend('Original','New','Inferred Primaries','Inferred Metastasis');
    %     elseif(n_infprim>0)
    %         legend('Original','New','Inferred Primaries');
    %     elseif(n_infmeta>0)
    %         legend('Original','New','Inferred Metastasis');
    %     else
    %         legend('Original','New');
    %     end
    %     title(sprintf('Day %i',data.t(i)));
    %     axis([1 3e7 0 2]);
end

% convert inferred distributions to same format as 'data' structure
% dist: [263×7 double]
% a: [263×1 double]
% t: [1 3 5 7 9 11 13]
% x: [263×1 double]
% cdf: [263×7 double]
% selector: {[48×1 double]  [57×1 double]  [50×1 double]  [35×1 double]  [24×1 double]  [19×1 double]  [18×1 double]}
% max: 2.8861e+07
% lowa: 0.4393
% lowx: 1.1713
% conv: [1×1 struct]
% xlabel: 'Tumor Size'
% ylabel: 'Number of Tumors per Body'
% plot_indices: [1 2 3 4 5 6 7]
% weights: [1 1 1 1 1 1 1]
xd = zeros(nxt,1);
% ad = zeros(nxt,1);
distd = zeros(nxt,nt);
cdfd = zeros(nxt,nt);
labelsd = zeros(nxt,nt);
selectord = cell(1,nt);
ninc = 0;
for i = 1:nt
    xd(ninc+1:ninc+length(temp_data(i).x)) = temp_data(i).x;
    distd(ninc+1:ninc+length(temp_data(i).x),i) = temp_data(i).dist;
    labelsd(ninc+1:ninc+length(temp_data(i).x),i) = temp_data(i).labels;
    ninc = ninc+length(temp_data(i).x);
end
% sort
[~,ord2] = sort(xd);
xd = xd(ord2);
ad = data.conv.v2a(xd);
distd = distd(ord2,:);
labelsd = labelsd(ord2,:);
% get selector and cdf
for i = 1:nt
    cdfd(:,i) = getStepIntegral(ad,distd(:,i));
    selectord{i} = find(distd(:,i)>0);
end
% place correct values into structure
inferred_data.dist = distd;
inferred_data.a = ad;
inferred_data.x = xd;
inferred_data.cdf = cdfd;
inferred_data.selector = selectord;
inferred_data.max = max(xd);
% store new format data into output
inferred_data.labels = labelsd;
inferred_data.alternate = temp_data;

end
