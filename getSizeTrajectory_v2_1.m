%% getSizeTrajectory_v2_1
%  Version 2.1
%  Author: Adeyinka Lesi
%  Date: 3/3/17
%  Project: Tumor Growth, Logarithmic Continuum Form


function [traj,t] = getSizeTrajectory_v2_1(key,x0s,dt,t0)
%  function [traj] = getSizeTrajectory_v1_0(key,x0)
%  key: struct of model details
%  x0s: vector, initial sizes
%  dt: scalar, time step
%  tf: scalar, final time
%  traj: 1xN vector, sizes
%  t: 1xN vector, times
%% VERSION HISTORY
%  2.1: x0 is now a vector

if(dt>0)
    t = 0:dt:t0;
else
    t = t0:dt:0;
end
traj = zeros(length(x0s),length(t));
traj(:,1) = x0s;
rates = key.RATES;

for i = 1:length(t)-1
    % fourth order runga-kutta to estimate mg_x
    xi = traj(:,i);
    dx1 = rates.growth(xi)-rates.death(xi)...
        -rates.shed(xi);
    est1 = max(1,xi+0.5*dt*dx1);
    dx2 = rates.growth(est1)-rates.death(est1)...
        -rates.shed(est1);
    est2 = max(1,xi+0.5*dt*dx2);
    dx3 = rates.growth(est2)-rates.death(est2)...
        -rates.shed(est2);
    est3 = max(1,xi+dt*dx3);
    dx4 = rates.growth(est3)-rates.death(est3)...
        -rates.shed(est3);
    traj(:,i+1) = max(1,xi+dt*(dx1+2*dx2+2*dx3+dx4)/6);
end
end

