function [avg_g,std_g] = getBackgroundMean_v1_0(g)
% remove outliers before recalculating mean
g_list = nonzeros(g);
avg1 = mean(g_list);
std1 = std(g_list);

g_list2 = g_list(g_list<avg1+3*std1);
avg_g = mean(g_list2);
std_g = std(g_list2);
end