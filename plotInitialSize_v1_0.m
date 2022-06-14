%% plotInitialSize_v1_0.m
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 10/15/20
%  Project: Fish Image Analysis

%% Version History
%  1.0: from plotInitialArea_v1_0.m - make it use fish data tables instead
%  this will be a plot used in paper

in = load('initial_area_zfall_032720.mat');
conv = getConverter();

imps = unique([in.silja_imps;in.zf2017_imps;in.zf2019_imps])';
meds = zeros(3,length(imps));

for i = 1:length(imps)
    meds(1,i) = median(in.silja_init_area1(in.silja_imps==imps(i)));
    meds(2,i) = median(in.zf2017_init_area1(in.zf2017_imps==imps(i)));
    meds(3,i) = median(in.zf2019_init_area1(in.zf2019_imps==imps(i)));
end

% figure;
% loglog(in.silja_imps,in.silja_init_area1,'o',in.zf2017_imps+0.25e5,...
%     in.zf2017_init_area1,'d',in.zf2019_imps+0.5e5,in.zf2019_init_area1,'s',...
%     'LineWidth',1.25)
% 
% legend({'Heilmann','Lesi2017','Lesi2019'},'Location','Best');
% xlabel('Inoculation Size (Cells)');
% ylabel('Day1 Total Area (px)');

figure;
loglog(in.silja_imps,conv.a2v(in.silja_init_area1),'o',in.zf2017_imps+0.25e5,...
    conv.a2v(in.zf2017_init_area1),'d',in.zf2019_imps+0.5e5,conv.a2v(in.zf2019_init_area1),'s',...
    'LineWidth',1.0);

xlabel('Inoculation Size (Cells)');
ylabel('Day 1 Estimated Size');
hold on;
% a1 = 10.^(1:0.1:5);
% loglog(areaToVol_v1_1(a1),a1);
% loglog(areaToVol_v2_0(a1),a1);
% loglog(imps+.75e5,meds,'*','MarkerSize',10);
hold off;
legend({'Heilmann','Lesi Batch 1','Lesi Batch 2',},'Location','Best');
axis([imps(1) imps(end)+0.5e5 -inf inf]);

% figure;
% loglog(in.silja_init_area1,in.silja_imps,'o',in.zf2017_init_area1,...
%     in.zf2017_imps+0.25e5,'d',in.zf2019_init_area1,in.zf2019_imps+0.5e5,'s',...
%     'LineWidth',1.0)
% 
% ylabel('Inoculation Size (Cells)');
% xlabel('Day1 Total Area (px)');
% hold on;
% a1 = 10.^(1:0.1:5);
% loglog(a1,areaToVol_v1_1(a1));
% loglog(a1,areaToVol_v2_0(a1));
% % loglog(imps+.75e5,meds,'*','MarkerSize',10);
% hold off;
% legend({'Heilmann','Lesi2017','Lesi2019','Old Estimate','New Estimate'},'Location','Best');
% axis([-inf inf imps(1) imps(end)+0.5e5]);
