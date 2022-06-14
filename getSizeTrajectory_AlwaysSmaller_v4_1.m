function [traj,t] = getSizeTrajectory_AlwaysSmaller_v4_1(key,x0,dt,t0)
%% getSizeTrajectory_v4_1
%  Version 4.1
%  Author: Adeyinka Lesi
%  Date: 11/7/20
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
%  3.1: adjustments to caculate values than are below 1
%  3.2: using log for differential equation
%  4.0: for time zero size calculator
%  4.1: changing how reduction is dealt with. Before, the code breaks out
%  of the loop if v(x)<0. Now, will simply set x(i+1)=x(i) when v(x(i))<0

if(dt>0)
    t = 0:dt:t0;
else
    t = linspace(t0,0,1+round(t0/abs(dt)));
end
ltraj = zeros(size(t));
ltraj(1) = log(x0);
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

for i = 1:length(ltraj)-1
    % check
    % fourth order runga-kutta to estimate mg_x
    dx1 = ((1-exp(ltraj(i))/cc)*rates.growth(exp(ltraj(i)))...
        -rfac(i)*rates.death(exp(ltraj(i)))...
        -rates.shed(exp(ltraj(i))))/exp(ltraj(i));
    est1 = max(lowx,ltraj(i)+0.5*dt*dx1);
    dx2 = ((1-exp(est1)/cc)*rates.growth(exp(est1))...
        -rfac(i)*rates.death(exp(est1))...
        -rates.shed(exp(est1)))/exp(est1);
    est2 = max(lowx,ltraj(i)+0.5*dt*dx2);
    dx3 = ((1-exp(est2)/cc)*rates.growth(exp(est2))...
        -rfac(i)*rates.death(exp(est2))...
        -rates.shed(exp(est2)))/exp(est2);
    est3 = max(lowx,ltraj(i)+dt*dx3);
    dx4 = ((1-exp(est3)/cc)*rates.growth(exp(est3))...
        -rfac(i)*rates.death(exp(est3))...
        -rates.shed(exp(est3)))/exp(est3);
    if((dx1+2*dx2+2*dx3+dx4) > 0 || dx1*dt < 0)
        % normal 
        ltraj(i+1) = max(lowx,ltraj(i)+dt*(dx1+2*dx2+2*dx3+dx4)/6);
    else
        ltraj(i+1) = ltraj(i);
    end
    %     elseif(abs(dt) <= low_tresh || abs(dx1*dt) <= low_tresh)
    %         % dt too low or
    %         % at balance point: no point in more calcs
    % %         fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e\n',dt,dx1,dx2,dx3,dx4);
    %         ltraj(i+1) = ltraj(i);
    %     else
    %         step_down_factor = 2;
    %         if(i > 1)
    %             % run with shorter time steps
    %             newt0 = t(i-1);
    %             newx0 = exp(ltraj(i-1));
    % %             fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e\n',dt,dx1,dx2,dx3,dx4);
    %             [newtraj,~] = getSizeTrajectory_AlwaysSmaller_v4_1(key,newx0,dt/step_down_factor,newt0);
    %
    %             ltraj(i:end) = log(newtraj((1:length(ltraj)-i+1)*step_down_factor+1));
    %         else
    %             % run with shorter time steps
    %             newt0 = t(i);
    %             newx0 = exp(ltraj(i));
    % %             fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e\n',dt,dx1,dx2,dx3,dx4);
    %             [newtraj,~] = getSizeTrajectory_AlwaysSmaller_v4_1(key,newx0,dt/step_down_factor,newt0);
    %
    %             ltraj(i+1:end) = log(newtraj((1:length(ltraj)-i)*step_down_factor+1));
    %         end
    %         break;
    %     end
end

traj = exp(ltraj);
end

