%% plotResults_v1_0
%  version 1.0
%  Author: Adeyinka Lesi
%  Date: 9/6/16
%  function [] = plotResults_v1_0(res)
%  res: struct, results from model and additional analysis

function [] = plotResults_v1_0(res)

key = res.key;
t = res.t;
tim = length(t);
pt = round(1+(tim-1)*0.25:(tim-1)*0.25:tim);
xp = res.xp;
xim = min(length(xp),res.max_index);

legendText = cell(1,length(pt));
for i = 1:length(pt)
    legendText{i} = ['t = ' num2str(t(pt(i))) ' Days'];
end

figure;
semilogy(xp(1:xim),res.dist(1:xim,pt));
title(key.TITLE);
ylabel('Population Density');
xlabel('Tumor Size');
legend(legendText,'Location','Best');

figure;
plot(xp(1:xim),res.cell_dist(1:xim,pt));
title(key.TITLE);
ylabel('Cell Density');
xlabel('Tumor Size');
legend(legendText,'Location','Best');

figure;
plot(t,res.num_tumors);
title(key.TITLE);
ylabel('Tumor Number')
xlabel('Time');

figure;
plot(t,res.num_meta);
title(key.TITLE);
ylabel('Metastasis Generatated')
xlabel('Time');