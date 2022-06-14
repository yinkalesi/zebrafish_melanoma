function key = getDataKey_v5_2(data,N_x,dt,name,pars,unk_dist)
%% getDataKey_v5_2
%  Version 5.2
%  Author: Adeyinka Lesi
%  Date: 1/25/21
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
%  4.0: using getKey_backtrack_v4_0
%  5.0: using getKey_backtrack_v5_0, which allows use of 
%  getTransformedDistribution_CTC_v1_0, which allows explicit CTC 
%  calculation and reseeding (requiring three new parameters)
%  5.1: using getTimeZeroSizes_v4_1
%  5.2: using STORED_RESULTS (default is to set to empty)

primary_init = max(data.x(data.selector{1}));
% primary_end = data.x(end);

% get key
key = getKey_backtrack_v5_2(pars,data.t(end),ceil(data.t(end)/dt),data,...
    N_x,primary_init,data.lowx,unk_dist,name);
key.TIMES_TO_STORE = data.t;
% key.HARD_TUMOR_SIZE_LIMIT = 10^(ceil(log10(primary_end))+1);

end