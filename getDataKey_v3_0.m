function key = getDataKey_v3_0(data,N_x,dt,name,pars,unk_dist)
%% getDataKey_v2_3
%  Version 2.3
%  Author: Adeyinka Lesi
%  Date: 10/09/17
%  function key = getDataKey_v2_2(data)
%  data: struct for experimental data
%  N_x: number of size intervals
%  dt: time increment
%  key: struct of model parameters
%% Version History
%  2.2: now adds new key field, key.TIME_ZERO_CALCULATOR
%  2.3: switched to getKey_backtrack_v2_0 - to implement carrying capacity
%  3.0: changes for intensity data, and changes to calculation of
%  time_zero_sizes (using getKey_backtrack_v3_0)

primary_init = max(data.x.*(data.dist(:,1)>0));
primary_end = data.x(end);

% get key
key = getKey_backtrack_v3_0(pars(1:2),pars(3:4),pars(5),pars(6:7),...
    data.t(1),data.t(end),ceil(data.t(end)/dt),data.x(data.selector{1}),...
    data.dist(data.selector{1},1),N_x,primary_init,data.lowx,unk_dist,name);
key.SOLUTION_FUNCTION = @getTransformedDistribution_timeSensitive_v1_0;
key.USING_METAGEN = 0;
key.TIMES_TO_STORE = data.t;
key.HARD_TUMOR_SIZE_LIMIT = 10^(ceil(log10(primary_end))+1);
key.CONV = data.conv;

end