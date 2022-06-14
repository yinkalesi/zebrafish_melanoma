function [x0,place_dist,place_rad,i_above1] = getTimeZeroSizes_v2_0(key,x_init,ts,lowx,unk_dist_key,unk_dist)
%% getTimeZeroSizes_v2_0
%  version 2.0
%  Author: Adeyinka Lesi
%  Date: 10/09/17

% get TIME_ZERO_SIZES using RK4 (step back in time using given velocities)
% need to handle tumors that are known to be present but are below
% detection limit. That implies that the specified size in x_init is
% irrelevant for these 'invisible' tumors and they will be placed based on
% a putative distribution

%% Version History
%  2.0: don't place tumors if growth time smaller than init t (needed to
%  switch to getSizeTrajectory_v3_0 for this)


if(isfield(key,'QUIET_MODE') && ~key.QUIET_MODE)
    display(['> ' mfilename]);
end
x0 = zeros(size(x_init));
place_dist = zeros(size(x_init));
place_rad = zeros(size(x_init));
if(ts>0)
    dt = ts/round(ts/key.TIME_STEP)/2;
    traj = getSizeTrajectory_v3_0(key,lowx,-dt,ts);
    lowx0 = traj(end);
    for j = 1:length(x_init)
        if(x_init(j) >= lowx)
            traj = getSizeTrajectory_v3_0(key,x_init(j),-dt,ts);
            x0(j) = max(0,traj(end));
        else
            x0(j) = 1+0.5*(lowx0-1);
            place_dist(j) = unk_dist_key.(unk_dist);
            place_rad(j) = 0.5*(lowx0-1);
        end
    end
else
    for j = 1:length(x_init)
        if(x_init(j) >= lowx)
            x0(j) = x_init(j);
        else
            x0(j) = 1+0.5*(lowx-1);
            place_dist(j) = unk_dist_key.(unk_dist);
            place_rad(j) = 0.5*(lowx-1);
        end
    end
end

% get rid of x0 values < 1
i_above1 = find(x0>=1);
if(isempty(i_above1))
    i_above1 = length(x0);
end
if(length(i_above1) < length(x0))
    fprintf('Using %i out of %i initial tumor values\n',length(i_above1),length(x0));
else
    fprintf('Using all %i initial tumor values\n',length(x0));
end
x0 = x0(i_above1);
place_dist = place_dist(i_above1);
place_rad = place_rad(i_above1);

