%% prePlotCalcs_v1_2
%  version 1.2
%  Author: Adeyinka Lesi
%  Date: 10/15/16
%  function [res] = prePlotCalcs_v1_2(res)
%  res: struct, results from model and additional analysis

function [res] = prePlotCalcs_v1_2(res)

xp = res.xp;
dist1 = res.dist;
dist2 = res.dist2;
dxs = res.key.FIELD.dxs;

cell_dist1 = zeros(size(dist1));

for i = 1:length(dist1(1,:))
    cell_dist1(:,i) = xp.*dist1(:,i)';
end

cell_dist2 = zeros(size(dist2));

for i = 1:length(dist2(1,:))
    cell_dist2(:,i) = xp.*dist2(:,i)';
end

dist_cum1 = zeros(size(dist1)+[1 0]);

for i = 1:length(dist1(1,:))
    % sum up the number of tumors cummulatively over tumor size
    for j = 1:length(dist1(:,1))
        dist_cum1(j+1,i) = dist1(j,i)*dxs(j)+ ...
                            dist_cum1(j,i);
    end
end

dist_cum2 = zeros(size(dist2)+[1 0]);

for i = 1:length(dist2(1,:))
    % sum up the number of tumors cummulatively over tumor size
    for j = 1:length(dist1(:,1))
        dist_cum2(j+1,i) = dist2(j,i)*dxs(j)+ ...
                            dist_cum2(j,i);
    end
end

cell_dist_cum1 = zeros(size(cell_dist1)+[1 0]);

for i = 1:length(cell_dist1(1,:))
    % sum up the number of cells cummulatively over tumor size
    for j = 1:length(cell_dist1(:,1))
        cell_dist_cum1(j+1,i) = cell_dist1(j,i)*dxs(j)+ ...
                                 cell_dist_cum1(j,i);
    end
end

cell_dist_cum2 = zeros(size(cell_dist2)+[1 0]);

for i = 1:length(cell_dist2(1,:))
    % sum up the number of cells cummulatively over tumor size
    for j = 1:length(cell_dist2(:,1))
        cell_dist_cum2(j+1,i) = cell_dist2(j,i)*dxs(j)+ ...
                                 cell_dist_cum2(j,i);
    end
end

res.cell_dist = cell_dist1;
res.dist_cum = dist_cum1;
res.cell_dist_cum = cell_dist_cum1;
res.cell_dist2 = cell_dist2;
res.dist_cum2 = dist_cum2;
res.cell_dist_cum2 = cell_dist_cum2;

%% Version History 

% 1.1: reflects changes to matrix sizes