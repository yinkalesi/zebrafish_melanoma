function [traj,t] = getSizeTrajectory_v3_0(key,x0,dt,t0)
%% getSizeTrajectory_v3_0
%  Version 3.0
%  Author: Adeyinka Lesi
%  Date: 7/21/17
%  Project: Tumor Growth, Logarithmic Continuum Form
%  function [traj] = getSizeTrajectory_v1_0(key,x0)
%  key: struct of model details
%  x0: scalar, initial size
%  dt: scalar, time step
%  tf: scalar, final time
%  traj: 1xN vector, sizes
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
if(dt>0)
    t = 0:dt:t0;
else
    t = t0:dt:0;
end
traj = zeros(size(t));
traj(1) = x0;
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
lowx = 1e-10; % minumum value for trajectory

for i = 1:length(traj)-1
    % check
    % fourth order runga-kutta to estimate mg_x
    dx1 = (1-traj(i)/cc)*rates.growth(traj(i))...
        -rfac(i)*rates.death(traj(i))...
        -rates.shed(traj(i));
    est1 = max(lowx,traj(i)+0.5*dt*dx1);
    dx2 = (1-est1/cc)*rates.growth(est1)...
        -rfac(i)*rates.death(est1)...
        -rates.shed(est1);
    est2 = max(lowx,traj(i)+0.5*dt*dx2);
    dx3 = (1-est2/cc)*rates.growth(est2)...
        -rfac(i)*rates.death(est2)...
        -rates.shed(est2);
    est3 = max(lowx,traj(i)+dt*dx3);
    dx4 = (1-est3/cc)*rates.growth(est3)...
        -rfac(i)*rates.death(est3)...
        -rates.shed(est3);
    traj(i+1) = max(lowx,traj(i)+dt*(dx1+2*dx2+2*dx3+dx4)/6);
end
end

