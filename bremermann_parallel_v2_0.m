function [x1,fval,flag,output] = bremermann_parallel_v2_0(func,x0,low,high,options)
%% bremermann_parallel_v2_0
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 04/26/19
%  Project: Parameter optimization
%  Description: algorithm picks a random direction and interpolates the 
%  function by 5-point lagrange interpolation in that direction to find 
%  minimums iteratively
%  Inputs:
%  func: function handle for function to be optimized
%  x0: 1xN vector of variables to be optimized
%  low*: 1xN vector of low bounds for variables
%  high*: 1xN vector of high bounds for variables
%  options*: stucture containing options for optimization
%           1) del: distance between interpolation points
%           2) maxIter: maximum number of iteraions (guesses)
%           3) TimeLimit: maximum time the algorithm can run (seconds)
%           4) Tolerance: end condition if function doesn't change by this
%           amount over a certain number of iterations
%           5) ToleranceLength: used to determine number of iterations
%           needed to test tolerance
%  Outputs:
%  res: 1xN vector of optimized variables
%  fval: scalar; minimum function value obtained
%  flag: integer; flag describing finishing state
%  output: struct; values describing run

%% Version History
%  1.0: enabling running func in parrallel using parfor
%  1.1: implementing changes from bremermann_v1_1
%  2.0: tracking tested parameters and fval vals to map parameter landscape

% set defaults
if(~exist('low','var'))
    low = -Inf*ones(1,length(x0));
end
if(~exist('high','var'))
    high = Inf*ones(1,length(x0));
end
if(~exist('options','var'))
    options = struct('del',0.05,'maxIter',1000,'TimeLimit',24*3600,...
        'Tolerance',1e-8,'ToleranceLength',10);
end

% del = 0.01 means each jump is 1% of the range of the parameter
range = high-low;
range(isinf(range)) = 1; % inf no limits supplied, default is 1

% timer
start = clock();

% set up parameters
lambda = [-2 -1 0 1 2]*options.del; % spacing for interpolation
N = length(lambda); % number of points
gap = range*options.del;
max_shift = floor(0.5*N);
endphrase = ['Finished after ' num2str(options.maxIter) ...
    ' iterations due to maxIter value'];
iter_count = options.maxIter;
flag = 0;
tolerance_count = 0;
tolerance_count_max = options.ToleranceLength*length(x0);

% calculate coefficients
precoef1 = zeros(1,N);
% explanation: these are denominators for Lagrange interpolation
for ip = 1:N
    precoef1(ip) = (-1)^(ip-1)*factorial(ip-1)*factorial(N-ip)*...
        (options.del)^(N-1);
end

precoef2 = zeros(N,N-1);
precoef2(:,1) = (N-1)*ones(N,1);
for row = 1:N
    % exclude current row
    sublam = lambda([1:row-1 row+1:end]);
    for col = 2:N-1
        % these values are needed to calculate interpolation coeficients
        % there are combinotorics involved as polynomials are being
        % multiplied - I am writing the code so that it works in the
        % arbitrary case (not just for N = 5)
        
        % explanation: to calculate precoef, need to add combinations
        % of values in sublam with no overlaps - the number of those
        % combinations is ncc (N-1 choose col-1). I am using each
        % unique ncc to find a unique combination, which I store in the
        % picks vector. Imagine a scenario where sublam = 1:4 and col-1
        % = 3; ncc = 4 and the combinations are 123,124,134,234 in
        % order of magnitude. The trick is to convert ni=2 into the
        % sequence 124 and so on...
        ncc = nchoosek(N-1,col-1);
        stepheight = N-1-(col-2);
        stepcounter = zeros(1,col-1);
        sum_jnoti = 0;
        for ni = 1:ncc
            picks = zeros(1,col-1);
            % make sure stepcounter is correct
            % 1) iterate lower counter if higher counter is full
            for ip = 1:length(picks)-1
                fpi = length(picks)-ip+1;
                if(stepcounter(fpi)>=stepheight)
                    stepcounter(fpi-1) = stepcounter(fpi-1) + 1;
                end
            end
            % 2) fix high counters that are full
            for ip = 1:length(picks)-1
                if(stepcounter(ip+1)>=stepheight)
                    stepcounter(ip+1) = stepcounter(ip);
                end
            end
            % get picks
            picks = (1:length(picks)) + stepcounter;
            % multiply and add to sum
            sum_jnoti = sum_jnoti + prod(sublam(picks));
            % iterate stepcounter
            stepcounter(end) = stepcounter(end) + 1;
        end
        precoef2(row,col) = (-1)^(col-1)*(N-col)*sum_jnoti;
    end
end
    

% iterate
x2 = x0;
x1 = x2;
fval1 = func(x2);
fval = fval1;
history = zeros(1,options.maxIter+1);
history(1) = fval1;
fout = zeros(1,options.maxIter*(2*N-2));
fout(1) = fval1;
fin = zeros(size(fout,2),length(x2));
fin(1,:) = x2;
neval = 1;

for iter = 1:options.maxIter
    % get random direction (normal distribution)
    % find limits to change:
    dlow = (x2-low)/gap;
    dhigh = (high-x2)/gap;
    % want to shift which x values are used when near edge
    shift = max(max_shift-floor(dlow),0)-max(max_shift-floor(dhigh),0);
    lambda2 = lambda + shift*options.del;
    % the range vectors is used as the std of the distribution
    dmax1 = abs((x2-low)/lambda2(1)/range);
    dmax2 = abs((high-x2)/lambda2(end)/range);
    dmax = min(dmax1,dmax2);
    dv = min(1,dmax1*0.5)*abs(randn(1,length(x2)));
    dv = min(dv,dmax);
    r = dv.*range;
    
    inters = zeros(1,N); % interpolation points
    inters_x = zeros(N,length(x2));
    % get inters_x
    for li = 1:N
        if(lambda2(li)~=0)
            % get x and make sure it is within the bounds (for the case of
            % non-uniform distribution)
            inters_x(li,:) = max(min(x2+lambda2(li)*r,high),low);
        else
            inters_x(li,:) = x2;
        end
    end
    
    % evaluate function
    parfor li = 1:N
        if(lambda2(li)~=0)
            inters(li) = func(inters_x(li,:));
        else
            inters(li) = fval1;
        end
    end
    Nlist = 1:N;
    Nlist = Nlist(lambda2~=0);
    nNlist = length(Nlist);
    fout(neval+(1:nNlist)) = inters(Nlist);
    fin(neval+(1:nNlist),:) = inters_x(Nlist,:);
    neval = neval+nNlist;
    
    % find minimum in polynomial interpolation (zeros of derivative)
    coefs = (inters./precoef1)*precoef2;
    if(any(isinf(coefs)) || any(isnan(coefs)))
        root_list = [];
    else
        % solve for zeros exactly (use root solver)
        root_list = roots(coefs) + shift*options.del;
    end
    % want real values only
    root_list = root_list(imag(root_list)==0);
    nrlist = length(root_list);
    fval_list = zeros(1,nrlist);
    x2_list = zeros(nrlist,length(x2));
    % find minimum by checking extrema
    
    % get list of extrema
    for rli = 1:nrlist
        x2_list(rli,:) = max(min(x2+root_list(rli)*r,high),low);
    end
    % evaluate function
    parfor rli = 1:nrlist
        fval_list(rli) = func(x2_list(rli,:));
    end
    fout(neval+(1:nrlist)) = fval_list;
    fin(neval+(1:nrlist),:) = x2_list;
    neval = neval+nrlist;
    
    % set fval1 to local internal minimum found (anything but current
    % value, unless current value is a root)
    all_fvals = [fval_list inters(lambda2~=0)];
    [~,min_i2] = min(all_fvals);
    fval1 = all_fvals(min_i2);
    all_x = [x2_list; inters_x(lambda2~=0,:)];
    x2 = all_x(min_i2,:);
    
    % find if there is a new global minimum
    if(fval >= fval1)
        % new minimum
        fval = fval1;
        x1 = x2;
    end
    
    % update history
%     disp([num2str(iter+1) ': ' num2str(fval1)]);
    history(iter+1) = fval1;
    % check tolerance
    if(abs(fval1-history(iter)) < options.Tolerance)
        tolerance_count = tolerance_count + 1;
        if(tolerance_count > tolerance_count_max)
            endphrase = ['Finished after ' num2str(iter) ...
                ' iterations due to Tolerance value'];
            flag = 1;
            iter_count = iter;
            break;
        end
    else
        % significant change: restart count
        tolerance_count = 0;
    end
    % check timer
    if(etime(clock(),start) > options.TimeLimit)
        endphrase = ['Finished after ' num2str(iter) ...
            ' iterations due to TimeLimit value'];
        flag = 0;
        iter_count = iter;
        break;
    end
end

output = struct('iterations',iter_count,'message',endphrase,...
    'totaltime',etime(clock(),start),'bounds',[low,high],...
    'history',history(1:iter_count+1),'inputs',fin(1:neval,:),...
    'outputs',fout(1:neval));