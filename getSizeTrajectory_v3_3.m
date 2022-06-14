function [traj,t] = getSizeTrajectory_v3_3(key,x0,dt,tval)
%% getSizeTrajectory_v3_0
%  Version 3.0
%  Author: Adeyinka Lesi
%  Date: 7/21/17
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
%  where t0 is just enough to get to size 1 and tval is too much; also
%  implementing death ramp
%  3.1: adjustments to caculate values than are below 1
%  3.2: using log for differential equation
%  3.3: calculate trajectories of multiple tumors at the same time; cc
%  factor uses sum of all tumors, removing step down option since this
%  version will not be used to get TIME_ZERO_SIZES

if(dt>0)
    t = 0:dt:tval;
else
    t = linspace(tval,0,1+round(tval/abs(dt)));
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
    dx1 = ((1-sum(exp(ltraj(:,i)))/cc).*rates.growth(exp(ltraj(:,i)))...
        -rfac(i)*rates.death(exp(ltraj(:,i)))...
        -rates.shed(exp(ltraj(:,i))))./exp(ltraj(:,i));
    est1 = max(lowx,ltraj(:,i)+0.5*dt*dx1);
    dx2 = ((1-sum(exp(est1))/cc).*rates.growth(exp(est1))...
        -rfac(i)*rates.death(exp(est1))...
        -rates.shed(exp(est1)))./exp(est1);
    est2 = max(lowx,ltraj(:,i)+0.5*dt*dx2);
    dx3 = ((1-sum(exp(est2))/cc).*rates.growth(exp(est2))...
        -rfac(i)*rates.death(exp(est2))...
        -rates.shed(exp(est2)))./exp(est2);
    est3 = max(lowx,ltraj(:,i)+dt*dx3);
    dx4 = ((1-sum(exp(est3))/cc).*rates.growth(exp(est3))...
        -rfac(i)*rates.death(exp(est3))...
        -rates.shed(exp(est3)))./exp(est3);
    
    % calculate next sizes
    ltraj(:,i+1) = max(lowx,ltraj(:,i)+dt*(dx1+2*dx2+2*dx3+dx4)/6);
    
%     if((dx1+2*dx2+2*dx3+dx4) > 0 || dx1 < 0)
%         % normal 
%     elseif(abs(dt) <= low_tresh || abs(dx1*dt) <= low_tresh)
%         % dt too low or
%         % at balance point: no point in more calcs
% %         fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e\n',dt,dx1,dx2,dx3,dx4);
%         ltraj(:,i+1:end) = ltraj(:,i)*ones(1,length(ltraj(:,i+1:end)));
%         break;
%     else
%         step_down_factor = 2;
%         if(i > 1)
%             % run with shorter time steps
%             newt0 = t(i-1);
%             newx0 = exp(ltraj(:,i-1));
% %             fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e\n',dt,dx1,dx2,dx3,dx4);
%             [newtraj,~] = getSizeTrajectory_v3_3(key,newx0,dt/step_down_factor,newt0);
%             
%             ltraj(:,i:end) = log(newtraj((1:length(ltraj)-i+1)*step_down_factor+1));
%         else
%             % run with shorter time steps
%             newt0 = t(i);
%             newx0 = exp(ltraj(:,i));
% %             fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e\n',dt,dx1,dx2,dx3,dx4);
%             [newtraj,~] = getSizeTrajectory_v3_3(key,newx0,dt/step_down_factor,newt0);
%             
%             ltraj(:,i+1:end) = log(newtraj((1:length(ltraj)-i)*step_down_factor+1));
%         end
%         break;
%     end
end

traj = exp(ltraj);
end

