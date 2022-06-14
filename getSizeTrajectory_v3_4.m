function [traj,t] = getSizeTrajectory_v3_4(key,x0,dt,t0,cap_factor)
%% getSizeTrajectory_v3_4
%  Version 3.4
%  Author: Adeyinka Lesi
%  Date: 7/3/18
%  Project: Tumor Growth, Logarithmic Continuum Form
%  function [traj] = getSizeTrajectory_v1_0(key,x0)
%  key: struct of model details
%  x0: 1xP vector, initial sizes
%  dt: scalar, time step
%  tf: scalar, final time
%  traj: PxN vector, sizes
%  t: 1xN vector, times
%% Version History
%  2.0: adjusted to work with negative dt with lowest possible size being 1
%  3.0: checks for carrying capacity and takes cc into account (assuming a
%  per fish basis for cc). Also, the lowest number traj can go down to is
%  now 1e-10 instead of 1 (want to avoid 0 because some of the rate
%  functions used will evaluate to NaN for 0 or negative numbers)
%  We want this lowest number to be less than one to distinguish the cases
%  where t0 is just enough to get to size 1 and t0 is too much; also
%  implementing death ramp
%  3.1: adjustments to caculate values than are below 1
%  3.2: using log for differential equation
%  3.3: calculate trajectories of multiple tumors at the same time; cc
%  factor uses sum of all tumors, removing step down option since this
%  version will not be used to get TIME_ZERO_SIZES
%  3.4: use the carrying capacity factor supplied (num_tumors/cc) - it is
%  assumed the time points match

if(dt>0)
    t = 0:dt:t0;
else
    t = linspace(t0,0,1+round(t0/abs(dt)));
end
ltraj = zeros(length(x0),length(t));
ltraj(:,1) = log(x0);
rates = key.RATES;

% check for carrying capacity
if(isfield(key,'CARRYING_CAPACITY'))
    cc = key.CARRYING_CAPACITY;
else
    cc = inf;
end
% check for ramping death rate
if(isfield(key,'USING_DEATH_RAMP'))
    rfac = key.DEATH_RAMP_FUNCTION(t);
else
    rfac = ones(1,length(t));
end
lowx = -inf; % minumum value for trajectory
low_tresh = 1e-4;

for i = 1:size(ltraj,2)-1
    % check
    % fourth order runga-kutta to estimate mg_x
    dx1 = ((1-cap_factor(i)).*rates.growth(exp(ltraj(:,i)))...
        -rfac(i)*rates.death(exp(ltraj(:,i)))...
        -rates.shed(exp(ltraj(:,i))))./exp(ltraj(:,i));
    est1 = max(lowx,ltraj(:,i)+0.5*dt*dx1);
    dx2 = ((1-0.5*(cap_factor(i)+cap_factor(i+1))).*rates.growth(exp(est1))...
        -rfac(i)*rates.death(exp(est1))...
        -rates.shed(exp(est1)))./exp(est1);
    est2 = max(lowx,ltraj(:,i)+0.5*dt*dx2);
    dx3 = ((1-0.5*(cap_factor(i)+cap_factor(i+1))).*rates.growth(exp(est2))...
        -rfac(i)*rates.death(exp(est2))...
        -rates.shed(exp(est2)))./exp(est2);
    est3 = max(lowx,ltraj(:,i)+dt*dx3);
    dx4 = ((1-cap_factor(i+1)).*rates.growth(exp(est3))...
        -rfac(i)*rates.death(exp(est3))...
        -rates.shed(exp(est3)))./exp(est3);
    
    % calculate next sizes
    ltraj(:,i+1) = max(lowx,ltraj(:,i)+dt*(dx1+2*dx2+2*dx3+dx4)/6);
    
end

traj = exp(ltraj);
end

