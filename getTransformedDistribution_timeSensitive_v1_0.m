function [result] = getTransformedDistribution_timeSensitive_v1_0(key)
%% getTransformedDistribution_timeSensitive_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 10/12/17
%  Project: Tumor Growth, Logarithmic Continuum Form
% getTransformedDistribution_v2_0 calculates the tumor size distribution in 
% an organism based on a tumor cell growth rate, death rate, and metastasis
% rate.  The formulas used are based on replacing the abscissa variable
% with its logorithm in order to facilitate extending the size of the
% domain.  
% params: a struct; contains an initial condition, domain, and time step
%   specifications for the calculations
% rates: a struct; contains rate functions
% result: a struct; contains the tumor size distributions, number of
%   tumors, and number of cells at each time
% function [result] = getTransformedDistribution_metagen_v2_0(params, rates)
%% Version History
%  1.0: generated from getTransformedDistribution_capacity_v1_1 - the
%  system matrix changes with time, but at a computational cost. The first
%  task is to include a ramping function for regression rate

params = key.FIELD;
rates = key.RATES;

if(~isfield(key,'QUIET_MODE'))
    quiet = 0;
else
    quiet = key.QUIET_MODE;
end

if(~quiet)
    disp(['> ' mfilename]);
end
start_time = clock(); %measure run time

% def: parameters for buffering feature
buffer = params.buffer;
reset_threshold = params.reset_threshold;

% define time variables
dt = params.dt;
% define some values used later
pm_init = params.init(1:end-1);
pm_limit = params.soft_limit;
pm_tmax = params.tmax;
if(isfield(key,'TIME_SHIFT'))
    tshift = key.TIME_SHIFT;
else
    tshift = 0;
end

if(params.death_start_time > 0)
    delay_death = params.using_death_starter;
else
    delay_death = 0;
end

% metagen implementation
if(isfield(key,'USING_METAGEN') && key.USING_METAGEN)
    % metastasis source will be used
    use_metagen = true;
    mg_start = key.METAGEN_START_SIZE;
    mg_end = key.METAGEN_END_SIZE;
    mg_tmax = key.METAGEN_MAX_TIME-tshift;
    if(delay_death)
        tmax = mg_tmax + params.death_start_time;
    else
        tmax = mg_tmax + pm_tmax;
    end
    
else
    % normal run
    use_metagen = false;
    if(delay_death)
        tmax = params.death_start_time-tshift;
    else
        tmax = pm_tmax-tshift;
    end
end
metagen_end_index = 0;

t = -tshift:dt:tmax;
% time points to store
if(~isfield(key,'SAVE_NUMBER'))
    key.SAVE_NUMBER = min(1000,length(t));
end
if(~isfield(key,'TIMES_TO_STORE'))
    key.TIMES_TO_STORE = tmax;
end
if(delay_death)
    dts1 = (tmax+tshift+pm_tmax-params.death_start_time)/(key.SAVE_NUMBER-1);
else
    dts1 = (tmax+tshift)/(key.SAVE_NUMBER-1);
end
ts1 = -tshift:dts1:tmax;
ts2 = key.TIMES_TO_STORE;
% identify indices for saved points
is1 = round((ts1+tshift)/dt)+1;
is2 = round((ts2+tshift)/dt)+1;
if(any(ts2+tshift > pm_tmax))
    warning('Cannot store all time points given');
    is2 = min(is2,length(t));
end

% how many time points in past to save
save_buffer = 5;

% define size variables
% dy = params.dy;
% dys = [params.dys params.dys(end)];
dxs = params.dxs;
y = params.y;
x = params.x;
xp = 0.5*(x(2:end)+x(1:end-1)); % value at grid centers
jmax = length(y);

% an appropriate correction to make the results better resemble the
% discrete case is to caculate the rate values a x-0.5 rather than x

% define rate functions
k_growth = rates.growth(x-0.5); % value at partition - 0.5
k_growth_p = rates.growth(xp-1.0); % shift of 1 is correction here

if(delay_death)
    k_death = zeros(1,length(x));
    k_death_p = zeros(1,length(x)-1);
    bc_rate_1 = 0;
    if(~quiet)
        disp('-> Death Starter will be used ...');
    end
else
    k_death = rates.death(x-0.5);
    k_death_p = rates.death(xp); % no shift here
    % death rate at boudary + correction terms (maybe)
    bc_rate_1 = -rates.death(x(1));
end
    

    
k_shed = rates.shed(x-0.5);
k_shed_p = rates.shed(xp); % no shift here

% k_meta_p = rates.meta(xp-0.5)+rates.meta_2deriv(xp-0.5)/24; % value at center - 0.5
k_meta_p = rates.meta(xp-0.5);

D = (k_growth_p + k_death_p + k_shed_p)/2; % value at function center
v = k_growth - k_death - k_shed; % value at partitions

% define carrying capacity

if(~isfield(key,'CARRYING_CAPACITY'))
    cap = inf;
else
    cap = key.CARRYING_CAPACITY;
end

% set up death rate ramp function
if(isfield(key,'USING_DEATH_RAMP') && key.USING_DEATH_RAMP)
    if(~quiet)
        disp('-> Death Ramp will be used ...');
    end
    c_ramp = 1-key.DEATH_RAMP_FUNCTION(t);
else
    c_ramp = zeros(1,length(t));
end
    
if(~quiet)
    disp('-> Creating Matrices ...');
end

% def: the following are coefficients for time discretization
trap_rule = dt*[-1/2 1/2];
implicit = dt*[-1 0];
explicit = dt*[0 1];

% chose schemes:
time_scheme = trap_rule;

% get matrices needed to solve system
[J,Jc,A,B,C] = getSystem(dxs,dt,time_scheme,jmax,bc_rate_1,D,v);

% construct  M vector
% the M vector contains parameters related to metastasis
M = zeros(1,jmax-1);
not1 = 2:jmax-1;
M(not1) = dt*k_meta_p(not1).*dxs(not1)/dxs(1); %added /dxs(1) for case when dxs(1)~=1

% construct integration templates (vectors used to integrate distribution)
tumorInt = dxs;
cellInt = dxs.*(xp-0.5);

% Process of obtaining solution at each time step:
% multiply B and p(t): O(n) since B is pentadiagonal (10n ops)
% multiply M and p(t) and add to B*p: O(n) (2n ops)
% get p(t+dt) using a tridiagonal matrix solver: O(n) (6n ops)

% p is a matrix containing distribution of tumor sizes (on a transformed scale)
% (row = tumor size) at different times (column = time)
p = zeros(jmax-1,save_buffer);
ps1 = zeros(jmax-1,length(ts1));
ps2 = zeros(jmax-1,length(ts2));

% metagen implementation
if(~use_metagen)
    % normal run
    % set the initial value
    p(:,1) = pm_init;
    soft_max = pm_limit;
    mg_size = 0;
else
    mg_size = zeros(length(mg_start),ceil(key.METAGEN_MAX_TIME/dt));
    mg_size(:,1) = mg_start;
    % change soft_max to increase speed
    soft_max = min(pm_limit,100);
end

if(~quiet)
    disp('-> Starting simulation ...');
end
% num_tumors contains the number of tumors at each time
num_tumors = zeros(1,length(t));
% This is calculated by trapezoid rule quadrature:
% The second term is a correction due to irregularity around initial size
num_tumors(1) = tumorInt*p(:,1);

% num_cells contains the number of cells at each time
num_cells = zeros(1,length(t));
% This is calculated by trapezoid rule quadrature:
num_cells(1) = cellInt*p(:,1);

% num_meta contains the number of metastasis at each time
num_meta = zeros(1,length(t));

% implementing metagen
if(use_metagen)
    if(mg_end-mg_size(:,1)>0)
        num_meta(1) = M*p(:,1)+dt*sum(rates.meta(mg_size(:,1)));
    else
        % gotten to mutagen end
        use_metagen = false;
        metagen_end_index = 1;
        % add in initial values
        p(:,1) = p(:,1)+pm_init;
        % need to change time array to end at pm_tmax
        t = (0:dt:pm_tmax)-tshift;
        if(~quiet)
            disp('-> Metagen end occured before any calculations');
        end
        % soft_max
        soft_max = max(soft_max,pm_limit);
        % generate metastasis
        num_meta(1) = M*p(:,1);        
    end
else
    % This will be calculated by multiplying the M vector and N(t)
    num_meta(1) = M*p(:,1);
end


% determine reporting frequency (progress will be printed as percentage)
if(isfield(key,'PERCENT_COMPLETE_INCREMENT'))
    repfreq = min((length(t)-1)*key.PERCENT_COMPLETE_INCREMENT,length(t)-1);
else
    repfreq = min(pm_tmax*0.1/dt,length(t)-2);
end
repfreq = max(1,repfreq);

fprintf('-> Progress: %7.2f%%', 0);

% vector contiaining soft_max value
max_loc = zeros(1,length(t));

% using while loop to enable loop to respond to length change of t
ms = 1:repfreq:(length(t)-repfreq);
mi = 1;
mi_max = length(ms);
ti_max = length(t)-1;

% next indices to be saved;
inis1 = 1;
inis2 = 1;
next_is1 = is1(inis1);
next_is2 = is2(inis2);

% record carrying capacity effects
cap_factor = zeros(1,length(t));
cap_factor(1) = num_cells(1)/cap;

while(mi <= mi_max)
        
    m = round(ms(mi));
    nmax = round(ms(mi)+repfreq)-1;
    mi = mi+1;
    
    for n = m:nmax
        % check for finish
        if(n > ti_max)
            if(~quiet)
                fprintf('\n');
                disp(['-> Stopping simulation early at t = ' num2str(t(ti_max)) ...
                    ' n = ' num2str(ti_max)]);
                fprintf('\n');
            end
            mi_max = mi-1;
            break;
        end
        % find index in compact matrix
        nc = mod(n-1,save_buffer)+1;
        ncp1 = mod(n,save_buffer)+1;
        
        % store data if neccesary
        if(n == next_is1)
            ps1(:,inis1) = p(:,nc);
            inis1 = inis1 + 1;
            if(inis1 <= length(is1))
                next_is1 = is1(inis1);
            end
        end
        % wait until after metagen to save this data
        if(~use_metagen)
            if(n == next_is2)
                ps2(:,inis2) = p(:,nc);
                inis2 = inis2 + 1;
                if(inis2 <= length(is2))
                    next_is2 = is2(inis2);
                end
            end
        end
       
        % reset soft_limit if values above threshold are within buffer
        while(p(soft_max,nc) > reset_threshold && soft_max < jmax-1)
            new_max = min(soft_max + buffer,jmax-1);
            if(new_max == jmax-1)
                if(~quiet)
                    fprintf('\n')
                    warning(['Buffer failure at t = ' num2str(t(n)) ...
                        ' n = ' num2str(n)]);
                    fprintf('-> Progress: %7.2f%%', 100*n/(length(t)-1));
                else
                    fprintf('\n')
                    disp(['Buffer failure at t = ' num2str(t(n)) ...
                          ' n = ' num2str(n)]);
                    fprintf('-> Progress: %7.2f%%', 100*n/(length(t)-1));
                end
            end
            
            % fprintf('\n-> soft_limit reset from %d to %d at t = %.2f, n = %d\n',soft_max,new_max,t(n),n);
            % fprintf('-> Progress: %7.2f%%', 100*n/(length(t)-1));
            soft_max = new_max;
        end
        max_loc(n) = soft_max;
        % multiply B (pentadiag matix) with N(t) and store as BN
        BCp = zeros(soft_max,1);
        % want to obtain better metastasis estimate
        % assuming we are using CN2 - we want num_meta(n+1/2)
        % will use extrapolation to approximate a value
        if(n > 1)
            nm_mid = (3*num_meta(n)-num_meta(n-1))/2; % using extrapolation for metastasis value...
        else
            nm_mid = num_meta(n);
        end  
        
%         D = (k_growth_p + k_death_p + k_shed_p)/2; % value at function center
%         v = k_growth - k_death - k_shed; % value at partitions
        % update system
        Dn = D-0.5*(k_growth_p*cap_factor(n)+k_death_p*c_ramp(n));
        vn = v-k_growth*cap_factor(n)+k_death*c_ramp(n);
        [J,Jc,A,B,C] = getSystem(dxs,dt,time_scheme,jmax,bc_rate_1,Dn,vn);
        
        % add metastasis to p(1,t+dt)
        Bpj = B(1,2:3)*p(1:2,nc);
        Cpj = C(1,2:3)*p(1:2,nc);
        BCp(1) = Bpj + Cpj + nm_mid; 
        for j = 2:(soft_max-1)
            Bpj = B(j,:)*p((j-1):(j+1),nc);
            Cpj = C(j,:)*p((j-1):(j+1),nc);
            BCp(j) = Bpj + Cpj;
        end
        % use no flux BC
        J_sm = J(soft_max,:);
        J_smp1 = zeros(size(J_sm));
        row = getMatrixRow(soft_max,dxs,J_sm,J_smp1,time_scheme);
        A_end = row(1:3);
        B_end = row(4:6);

        BCp(soft_max) = B_end(1:2)*p((soft_max-1):soft_max,nc);
        
        % solve using tridiagonal matrix algorithm
        p(1:soft_max,ncp1) = solveTridiag([A(1:soft_max-1,:); A_end], BCp);
        
        % update other variables
        num_tumors(n+1) = tumorInt(1:soft_max)*p(1:soft_max,ncp1);
        num_cells(n+1) = cellInt(1:soft_max)*p(1:soft_max,ncp1);        
        cap_factor(n+1) = num_cells(n+1)/cap;
        
        % check that metagen max time has not been passed
        if(use_metagen && t(n+1)>mg_tmax) 
            % stop simulation
%             use_metagen = false;
%             metagen_end_index = n;
%             p(:,nc) = p(:,nc)+pm_init;
%             p(:,ncp1) = p(:,ncp1)+pm_init;
%             change time array
%             tshift_mg = t(metagen_end_index);
%             ipm_tmax = ceil((pm_tmax+tshift+tshift_mg)/dt)+1;
%             t = t(1:ipm_tmax)-tshift_mg;
%             t = [t(1:metagen_end_index-1)-tshift_mg 0:dt:pm_tmax+dt];
%             if(t(end-1) >= pm_tmax)
%                 t = t(1:end-1);
%             end
%             ti_max = length(t)-1;
%             ts1 = ts1-tshift_mg;
%             update TIME_TO_SAVE points
%             is2 = round((ts2+tshift+tshift_mg)/dt)+1;
%             next_is2 = is2(inis2);
%             save if necessary
%             if(n == next_is2)
%                 ps2(:,inis2) = p(:,nc);
%                 inis2 = inis2 + 1;
%                 if(inis2 <= length(is2))
%                     next_is2 = is2(inis2);
%                 end
%             end
%             
%             if(~quiet)
%                 fprintf('\n');
%                 disp('-> Time shifted due to metagen end');
%                 fprintf(['Time shift at t = ' num2str(-t(1)-tshift) ...
%                     ' n = ' num2str(n) '\n']);
%                 fprintf('-> Progress: %7.2f%%', 100*n/(length(n)-1));
%             end
%             soft_max = max(soft_max,pm_limit);
            fprintf('\n');
            warning('Metagen Time Limit Reached - Simulation will cease');
            fprintf('-> Progress: %7.2f%%', 100*n/(length(t)-1));
            ti_max = n;
            break;
        end
        
        % implementing metagen
        if(use_metagen)
            % fourth order runga-kutta to estimate mg_x
            dmg_x1 = rates.growth(mg_size(:,n))-rates.death(mg_size(:,n))...
                -rates.shed(mg_size(:,n));
            mg_est1 = mg_size(:,n)+0.5*dt*dmg_x1;
            dmg_x2 = rates.growth(mg_est1)-rates.death(mg_est1)...
                -rates.shed(mg_est1);
            mg_est2 = mg_size(:,n)+0.5*dt*dmg_x2;
            dmg_x3 = rates.growth(mg_est2)-rates.death(mg_est2)...
                -rates.shed(mg_est2);
            mg_est3 = mg_size(:,n)+dt*dmg_x3;
            dmg_x4 = rates.growth(mg_est3)-rates.death(mg_est3)...
                -rates.shed(mg_est3);
            mg_size(:,n+1) = mg_size(:,n)+dt/6*(dmg_x1+2*dmg_x2+2*dmg_x3+dmg_x4);
            if(mg_size(:,n+1) < 0)
                fprintf('\n');
                warning('Inappropriate parameter set - Simulation will cease');
                fprintf('-> Progress: %7.2f%%', 100*n/(length(t)-1));
                ti_max = n;
                break;
            end
            if(mg_end-mg_size(:,n+1)>0)
                num_meta(n+1) = M(1:soft_max)*p(1:soft_max,ncp1) ...
                    + dt*sum(rates.meta(mg_size(:,n+1)));
            else
                % gotten to mutagen end
                use_metagen = false;
                metagen_end_index = n;
                % add in initial values
                p(:,nc) = p(:,nc)+pm_init;
                p(:,ncp1) = p(:,ncp1)+pm_init;
                % change time array
                tshift_mg = t(metagen_end_index);
                ipm_tmax = ceil((pm_tmax+tshift+tshift_mg)/dt)+1;
                t = t(1:ipm_tmax)-tshift_mg;
%                 t = [t(1:metagen_end_index-1)-tshift_mg 0:dt:pm_tmax+dt];
%                 if(t(end-1) >= pm_tmax)
%                     t = t(1:end-1);
%                 end
                ti_max = length(t)-1;
                ts1 = ts1-tshift_mg;
                % update TIME_TO_SAVE points
                is2 = round((ts2+tshift+tshift_mg)/dt)+1;
                next_is2 = is2(inis2);
                % save if necessary
                if(n == next_is2)
                    ps2(:,inis2) = p(:,nc);
                    inis2 = inis2 + 1;
                    if(inis2 <= length(is2))
                        next_is2 = is2(inis2);
                    end
                end
            
                if(~quiet)
                    fprintf('\n');
                    disp('-> Time shifted due to metagen end');
                    fprintf(['Time shift at t = ' num2str(-t(1)-tshift) ...
                        ' n = ' num2str(n) '\n']);
                    fprintf('-> Progress: %7.2f%%', 100*n/(length(t)-1));
                end
                % soft_max
                soft_max = max(soft_max,pm_limit);
                % generate metastasis
                num_meta(n+1) = M(1:soft_max)*p(1:soft_max,ncp1);
            end
        else
            num_meta(n+1) = M(1:soft_max)*p(1:soft_max,ncp1);
        end
        
    end
    fprintf([repmat(8,1,8) '%7.2f%%'], 100*nmax/(length(t)-1));
end

% store last time point if needed
n = length(t);
nc = mod(n-1,save_buffer)+1;

% store data if neccesary
if(n == next_is1)
    ps1(:,inis1) = p(:,nc);
    inis1 = inis1 + 1;
end
% wait until after metagen to save this data
if(~use_metagen)
    if(n == next_is2)
        ps2(:,inis2) = p(:,nc);
        inis2 = inis2 + 1;
    end
end
 
fprintf('\n');
if(~quiet)
    disp('-> End of Simulation');
end

run_time = etime(clock(),start_time);

t = t(1:is1(inis1-1));
ts1 = ts1(1:inis1-1);
ts2 = ts2(1:inis2-1);
ps1 = ps1(:,1:inis1-1);
ps2 = ps2(:,1:inis2-1);
num_tumors = num_tumors(1:is1(inis1-1));
num_cells = num_cells(1:is1(inis1-1));
num_meta = num_meta(1:is1(inis1-1));
max_loc = max_loc(1:is1(inis1-1));
cap_factor = cap_factor(1:is1(inis1-1));
c_ramp = c_ramp(1:is1(inis1-1));

% check if using death starter
if(delay_death && params.death_start_time < pm_tmax)
    % check for special case where surgery is also performed
    if(params.using_surgery && params.surg_t == params.death_start_time)
        % do nothing, let the simulation restart
    else 
        newpams = params;
        newpams.using_death_starter = false;
        newpams.death_start_time = -1;
        newpams.tmax = pm_tmax-params.death_start_time;
        newpams.init = [ps1(:,inis1-1); 0];
        newpams.soft_limit = soft_max;
        
        newkey = key;
        newkey.FIELD = newpams;
        newkey.TIME_SHIFT = -params.death_start_time;
        newkey.TIMES_TO_STORE = key.TIMES_TO_STORE(key.TIMES_TO_STORE>...
            params.death_start_time);%-params.death_start_time;
        newkey.USING_METAGEN = false;
        newkey.SAVE_NUMBER = key.SAVE_NUMBER - inis1 + 2;
        
%         continue running simulation to obtain final results
        thisFunc = str2func(mfilename);
        finres = thisFunc(newkey);
        
        if(ts1(inis1-1)==finres.t(1))%+params.death_start_time)
            ts1 = [ts1(1:inis1-2) (finres.t)];%+params.death_start_time)];
            ps1 = [ps1(:,1:inis1-2) finres.dist];
            t = [t(1:end-1) (finres.t3)];%+params.death_start_time)];
            num_tumors = [num_tumors(1:end-1) finres.num_tumors];
            num_cells = [num_cells(1:end-1) finres.num_cells];
            num_meta = [num_meta(1:end-1) finres.num_meta];
            max_loc = [max_loc(1:end-1) finres.max_loc];
            cap_factor = [cap_factor(1:end-1) finres.capacity_factor];
            c_ramp = [c_ramp(1:end-1) finres.complementary_ramp];
        else
            ts1 = [ts1(1:inis1-1) (finres.t)];%+params.death_start_time)];
            ps1 = [ps1(:,1:inis1-1) finres.dist];
            t = [t (finres.t3)];%+params.death_start_time)];
            num_tumors = [num_tumors finres.num_tumors];
            num_cells = [num_cells finres.num_cells];
            num_meta = [num_meta finres.num_meta];
            max_loc = [max_loc finres.max_loc];
            cap_factor = [cap_factor finres.capacity_factor];
            c_ramp = [c_ramp finres.complementary_ramp];
        end
        if(ts2(inis2-1)==finres.t2(1))%+params.death_start_time)
            ts2 = [ts2(1:inis2-2) (finres.t2)];%+params.death_start_time)];
            ps2 = [ps2(:,1:inis2-2) finres.dist2];
        else
            ts2 = [ts2(1:inis2-1) (finres.t2)];%+params.death_start_time)];
            ps2 = [ps2(:,1:inis2-1) finres.dist2];
        end
        run_time = run_time + finres.run_time;
        soft_max = finres.max_index;
        k_death = finres.death_rate;
    end
end
    
% create results struct
result = struct('t',ts1,'dist',ps1,...
    't2',ts2,'dist2',ps2,...
    't3',t,'num_tumors',num_tumors,...
    'num_cells',num_cells,...
    'num_meta',num_meta, ...
    'run_time',run_time, 'y', y, 'M',M, 'xp', xp, ...
    'x',x,'growth_rate',k_growth,'death_rate',k_death,...
    'shed_rate',k_shed,'meta_rate',k_meta_p,'max_index',soft_max,...
    'max_loc',max_loc,'metagen_size',mg_size,...
    'metagen_end_index',metagen_end_index,...
    'capacity_factor',cap_factor,'complementary_ramp',c_ramp);

if(~quiet)
    disp(['-> Run Time: ' num2str(result.run_time) ' s']);
end

end

function row = getMatrixRow(i,dxs,J1,J2,tcoef)
fluxDiff = [0 J2] - [J1 0];

arow = [0 1 0] - tcoef(1)/dxs(i)*fluxDiff(1:3);
brow = [0 1 0] - tcoef(2)/dxs(i)*fluxDiff(1:3);

row = [arow brow];
end

function corr = getMatrixCorr(i,dxs,Jc1,Jc2,dt)
    corr = dt*([Jc1 0]-[0 Jc2])/dxs(i);
end

function [flux] = getFlux(i,dxs,D,v)
% Ji = v*pi - d(D*pi)/dx

denom = 2/(dxs(i-1)+dxs(i));
flux = zeros(1,2);
    
if(v(i) < 0)
    flux(2) = v(i)-D(i)*denom;
    flux(1) = D(i-1)*denom;
elseif(v(i) > 0)
    flux(2) = -D(i)*denom;
    flux(1) = v(i)+D(i-1)*denom;
else
%   don't need velocity terms here since v(i)=0
    flux(2) = -D(i)*denom;
    flux(1) = D(i-1)*denom;
end
end

function [corr] = getFluxCorr(i,dxs,v)
denom = 1/(dxs(i-1)+dxs(i));
if(v(i) > 0)
    corr = v(i)*denom*[-1 1];
else
    corr = v(i)*denom*[1 -1];
end
end

function [J,Jc,A,B,C] = getSystem(dxs,dt,time_scheme,jmax,bc_rate_1,D,v)
% System is: A*N(n+1) = (B+M)*N(n) where n is the time index
% construct A and B matrices
% A is tridiagonal, so only need three columns
A = zeros(jmax-1,3);
% B is pentadiagonal, so five columns
B = zeros(jmax-1,3);
% C is a correction matrix (2 columns)
C = zeros(jmax-1,3);

% def: J is a matrix containing the discretization values for the flux
J = zeros(jmax,2);
% def: Jc is for flux correction
Jc = zeros(jmax,2);

% flux out of boundary (tumor removal)
% getFlux(i,gaps,growth,k_death,k_shed)
J(1,:) = [0 bc_rate_1]; % removed /dxs(1)
J(2,:) = getFlux(2,dxs,D,v);
Jc(1,:) = zeros(size(Jc(1,:)));
Jc(2,:) = getFluxCorr(2,dxs,v);

% getMatrixRow(i,dxs,J1,J2,tcoef)
row = getMatrixRow(1,dxs,J(1,:),J(2,:),time_scheme);
A(1,:) = row(1:3);
B(1,:) = row(4:6);
C(1,:) = getMatrixCorr(1,dxs,Jc(1,:),Jc(2,:),dt);

% fill in second row
J(3,:) = getFlux(3,dxs,D,v);
Jc(3,:) = getFluxCorr(3,dxs,v);

% getMatrixRow(i,dxs,J1,J2,tcoef)
row = getMatrixRow(2,dxs,J(2,:),J(3,:),time_scheme);
A(2,:) = row(1:3);
B(2,:) = row(4:6);
C(2,:) = getMatrixCorr(2,dxs,Jc(2,:),Jc(3,:),dt);

% fill in other rows until jmax-3
for i = 3:(jmax-3)
    % fill in ith row
    J(i+1,:) = getFlux(i+1,dxs,D,v);
    Jc(i+1,:) = getFluxCorr(i+1,dxs,v);

    % getMatrixRow(i,dxs,J1,J2,tcoef)
    row = getMatrixRow(i,dxs,J(i,:),J(i+1,:),time_scheme);
    A(i,:) = row(1:3);
    B(i,:) = row(4:6);
    C(i,:) = getMatrixCorr(i,dxs,Jc(i,:),Jc(i+1,:),dt);
end

% fill in jmax-2 row (since near edge, forcing use of fd1 scheme)
i = jmax-2;
J(i+1,:) = getFlux(i+1,dxs,D,v);
Jc(i+1,:) = getFluxCorr(i+1,dxs,v);

% getMatrixRow(i,dxs,J1,J2,tcoef)
row = getMatrixRow(i,dxs,J(i,:),J(i+1,:),time_scheme);
A(i,:) = row(1:3);
B(i,:) = row(4:6);
C(i,:) = getMatrixCorr(i,dxs,Jc(i,:),Jc(i+1,:),dt);

% fill in last row 
% BC options:
% 1) No flux: J(jmax+1) = 0 *(in use)
% 2) Pass through: J(jmax) = J(jmax+1)
i = jmax-1;
J(i+1,:) = zeros(size(J(i,:)));
Jc(i+1,:) = zeros(size(Jc(i,:)));

% getMatrixRow(i,dxs,J1,J2,tcoef)
row = getMatrixRow(i,dxs,J(i,:),J(i+1,:),time_scheme);
A(i,:) = row(1:3);
B(i,:) = row(4:6);
C(i,:) = getMatrixCorr(i,dxs,Jc(i,:),Jc(i+1,:),dt);

end
%% Version Log
%  1.0: named code getTransformedDistribution_metagen_v2_0. Added
%  capability to generate metastasis from a source not on the grid (to
%  allow metastasis generation following a known tragectory)
%  1.1: streamlining the code to make more efficient; boundary condition is
%  now no flux;
%  2.0: only storing a select number of variables and adding specific
%  storage a subset of timepoints to facilitate data comparison

% Logs for getTransformedDistribution
% 1.0: Speed improvement through soft_limit, using predefined rate
% functions and derivatives to spatial discretization
% [RK2+CN/(FD2,BD2,CD)]
% 1.1: Removed use of trapezoid rule; might be causing accuracy issues.
% There was an error with redundant k0 value in BC that was addressed
% [RK2/(FD2,BD2,CD)]
% 1.2: Using only forward difference, first order; found error: the code is
% actually calculating N(:,t+3*dt/2) rather than N(:,t+dt)
% [RK2/(FD1)]
% 1.3: Not running Runge-Kutta
% [CN/(FD2,BD2,CD)]
% 2.0: Implementing irregular mesh around initial tumor size so that can
% have the same starting point for irregular tumor sizes
% Noticed and fixed error with calculating M matrix (added dy)
% (7/27/15: reorganized code to simplify)
% (7/29/15: fixed an error in irregular mesh discretization scheme)
% 3.0: Enabling multiple initial tumors
% 4.0: Enabling generic descretization scheme
% (8/14/15): Found an error in tridiagonal matrix solver that caused issues
% with the boundary condition
% [(CN,EE)/(FD2,BD2,CD)]; BC:[(CN,IE)/(FD1)]
% 4.1: Different discretization scheme
% [(CN,EE)/(FD2,BD2,CD)]; BC:[(CN)/(FD1)]
% 4.2: Different discretization scheme
% [(CN,EE)/(FD2,BD2,CD)]; BC:[(CN,EE)/(FD1)]
% 4.3: Using first order (but hopefully more stable) scheme
% [(CN)/(FD1,BD1,CD)]; BC:[(CN)/(FD1)]
% 5.0: Adjusting formula for calculating number of tumors and number of
% cells (numberical integration)
% [(CN)/(FD1,BD1,CD)]; BC:[(CN)/(FD1)]
% 5.1: A correction needs to be made to boundary conditions so it gives the
% proper value when dy is not 1.  The zero order derivative part at the
% boundary should be adjusted
% [(CN)/(FD1,BD1,CD)]; BC:[(CN)/(FD1)]
% 6.0: Boundary condition overhaul: N1 is treated as a discrete part of the
% system and it is continuous for x>2
% (8/24/16): using dx instead of dy to calculate integrals
% 6.1: RK2
% 6.3: trying to eliminate apparent error
% 7.0: Implementation of finite volume style scheme to conserve tumor
% number (bad velocity accuracy)
% 7.1: Making some corrections to hopefully improve algorithm;
% however, this one has bad diffusion accuracy
% 7.2: Another go at an improvement - the idea is to use finite volume to
% solve the continuous problem rather than trying to use it, or at least
% the idea of fluxes, to solve the discrete problem (done for 7.1)
% in addition defined a velocity and a diffusion coeficient instead of
% using the growth and death rates explicitly
% 7.3: Fixing an issue where its is possible for d(pi)/dt < 0 while pi = 0
% 7.5: Trying to fix stability issue
% 7.6: Error in formula: D = (k1+f)/2.
% 8.0: Error: Metastasis vector must be multiplied by dt
% 8.1: Changing discretization expression for fluxes; hopefully will give
% more stable simulations; implemented upwind scheme where only upstream
% terms are used in convection. 
% 8.2: Implemented DEATH_START_TIME
% 8.3: improving metastasis accuracy - first, let's try extrapolation and
% then use simpson's rule integration
% 8.4: possible error in handling metastasis that shows up if dxs(1)~=1
% 9.0: 1) partitions at integer points
%      2) function values on half-intervals
%      3) initial distribution as a Gaussian curve or 
%      4) input is now key
% 9.1: Adding additional term (derivative of v(x)) for more accuracy. Now
% all second order terms are included
% 9.2: Explicit version, for testing
% 9.3: As many correction terms as possible
%      1) metastasis term corrections
%      2) shift rate values
%      The use of upwind creates an error
% 9.4: not using upwind (much higher accuracy)
% 9.5: using upwind with correction terms to compensate (the corrections
% are added in an explicit fashion which seems to improve stabililt without
% damaging stability (for dx >= 1, results for dx < 1 seem to start going
% to infinity) - It turns out that this version does produce negative
% values; its just that the negative values don't tend to blow up
% 9.6: modifying the correction to prevent negative values from appearing