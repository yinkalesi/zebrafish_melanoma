%% paper1_PlotCompilations_v1_0
%  Version 1.4
%  Author: Adeyinka Lesi
%  Contributors: Dr. Silja Heilmann
%  Date: 1/6/22
%  Project: Experimental Paper
%% Version History
%  1.4: changed range of data

startt = clock();

savedir = 'C:\Users\Yinka\Documents\CCNY\Rumschitzki Group\Continuous Tumor Growth Model\paper1_plots_011622\';

%% 1) plot of fish data distribution with optimal parameters
% (make sure to mark inferred data)
mrs_dat=load('ct10_lt015_ht04_021720_rev2_fits_083120_highlow_runModelHelper_zf_fits_v2_1.mat');
fontsizes = [6 8 7 6]*2;

mrs_dat.mrs{1,1}.fitter.objective(mrs_dat.mrs{1,1},mrs_dat.mrs{1,1}.data);
axis([10 3e7 0 max(mrs_dat.mrs{1,1}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('');
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sf_r_cdf',savedir),'-djpeg');

plotMetaOriginTimes_v1_1(mrs_dat.mrs{3,1},mrs_dat.mrs{1,1},'Male','Female');
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Days Post Inoculation','FontSize',2*fontsizes(2));
ylabel('Number of Metastasis','FontSize',2*fontsizes(3));
title('');
legend('Location','NorthWest','FontSize',2*fontsizes(4));
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',2*fontsizes(1));
print(sprintf('%smf_r_met',savedir),'-djpeg');

mrs_dat.mrs{3,1}.fitter.objective(mrs_dat.mrs{3,1},mrs_dat.mrs{3,1}.data);
axis([10 3e7 0 max(mrs_dat.mrs{3,1}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sm_r_cdf',savedir),'-djpeg');

mrs_dat.mrs{2,1}.fitter.objective(mrs_dat.mrs{2,1},mrs_dat.mrs{2,1}.data);
axis([10 3e7 0 max(mrs_dat.mrs{2,1}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sf_nr_cdf',savedir),'-djpeg');

mrs_dat.mrs{4,1}.fitter.objective(mrs_dat.mrs{4,1},mrs_dat.mrs{4,1}.data);
axis([10 3e7 0 max(mrs_dat.mrs{4,1}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sm_nr_cdf',savedir),'-djpeg');

% Michaelis-Menten immComp plots
mm_mrs_dat=load('ct10_lt015_ht04_021720_final_mm_fits_123120_runModelHelper_mm_v1_2.mat');

mm_mrs_dat.mrs{2,3}.fitter.objective(mm_mrs_dat.mrs{2,3},mm_mrs_dat.mrs{2,3}.data);
axis([10 3e7 0 max(mm_mrs_dat.mrs{2,3}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sf_mm_cdf',savedir),'-djpeg');

mm_mrs_dat.mrs{4,2}.fitter.objective(mm_mrs_dat.mrs{4,2},mm_mrs_dat.mrs{4,2}.data);
axis([10 3e7 0 max(mm_mrs_dat.mrs{4,2}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sm_mm_cdf',savedir),'-djpeg');

% 5e6 inoculation long term plots
mm_mrs_dat=load('ct10_lt015_ht04_021720_mm_FisZ_inoc5e6V_0040821_fit1_runModelHelper_zf_fits_v2_9.mat');

mm_mrs_dat.mrs{5,1}.fitter.objective(mm_mrs_dat.mrs{5,1},mm_mrs_dat.mrs{5,1}.data);
axis([10 3e7 0 max(mm_mrs_dat.mrs{5,1}.data.cdf(end,:))]);
[~,~,zf19leg_h,~] = legend('Location','BestOutside','FontSize',fontsizes(4));
zf19leg_sel_beg = [1:3 5:6 8];
zf19_leg_sel_rest = 9:length(zf19leg_h);
zf19_leg_sel_rest = zf19_leg_sel_rest(mod(zf19_leg_sel_rest,4)==1|mod(zf19_leg_sel_rest,4)==2);
zf19leg_sel = [zf19leg_sel_beg zf19_leg_sel_rest];
zf19data = mm_mrs_dat.mrs{5,1}.data;
zf19leg_text = cell(1,4*length(zf19data.plot_indices));
for i = 1:length(zf19data.plot_indices)
    zf19leg_text{4*i-3} = ['Data, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
    zf19leg_text{4*i-2} = ['Model, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
    zf19leg_text{4*i-1} = ['Inf. Prim, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
    zf19leg_text{4*i-0} = ['Inf. Meta, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
end
legend(zf19leg_h(zf19leg_sel),zf19leg_text(zf19leg_sel),'Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sf_mm_inoc5e6_cdf',savedir),'-djpeg');

mm_mrs_dat.mrs{6,1}.fitter.objective(mm_mrs_dat.mrs{6,1},mm_mrs_dat.mrs{6,1}.data);
axis([10 3e7 0 max(mm_mrs_dat.mrs{6,1}.data.cdf(end,:))]);
[~,~,zf19leg_h,~] = legend('FontSize',4,'Location','Best');
zf19leg_sel_beg = [1:3 5:6 8];
zf19_leg_sel_rest = 9:length(zf19leg_h);
zf19_leg_sel_rest = zf19_leg_sel_rest(mod(zf19_leg_sel_rest,4)==1|mod(zf19_leg_sel_rest,4)==2);
zf19leg_sel = [zf19leg_sel_beg zf19_leg_sel_rest];
zf19data = mm_mrs_dat.mrs{6,1}.data;
zf19leg_text = cell(1,4*length(zf19data.plot_indices));
for i = 1:length(zf19data.plot_indices)
    zf19leg_text{4*i-3} = ['Data, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
    zf19leg_text{4*i-2} = ['Model, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
    zf19leg_text{4*i-1} = ['Inf. Prim, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
    zf19leg_text{4*i-0} = ['Inf. Meta, t = ' num2str(1+zf19data.t(zf19data.plot_indices(i)))];
end
legend(zf19leg_h(zf19leg_sel),zf19leg_text(zf19leg_sel),'Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('')
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%sm_mm_inoc5e6_cdf',savedir),'-djpeg');

%% plot of Iwata distribution
iwata_res1 = load('iwata_results_P1_010621_runModelHelper_iwata_vP1_save.mat');
iwata_res2 = load('iwata_results_P2_010621_runModelHelper_iwata_vP2_save.mat');

hpic = 3.125*61*2;
wpic = 3.125*61*2;
legpos4 = [0.6 0.7 0.2 0.2];
xtic1 = 10.^(7:1:12);
tic_def = 'auto';

plotCDF2_vP1(iwata_res1.mg_resf3,iwata_res2.mg_resf3b,iwata_res1.dataf3,4:6,[1e7 3e11 0 50],round(hpic),round(wpic),'Tumor Size','Cumulative Number of Tumors',fontsizes,legpos4,xtic1,tic_def);
print(sprintf('%siwata_cdf',savedir),'-djpeg');

plotCellDist2_vP1(iwata_res1.mg_resf3,iwata_res2.mg_resf3b,iwata_res1.dataf3,4:6,[1e7 3e11 0 4],hpic,wpic,'Tumor Size','Size-Weighted Number Density',fontsizes,legpos4,tic_def,tic_def);
print(sprintf('%siwata_wpdf',savedir),'-djpeg');

%% plot of Silja data distribution (need refit for metastasis portion)
silja_dat=load('Silja_022520_r2_s1_032221_ap_040421_bestfit2_runModelHelper_Silja_v1_4.mat');
silja_gomp = load('Silja_022520_r2_s1_031020_ap_GompSil_021921_runModelHelper_GompSilja_v1_0.mat');

silja_dat.mrs{1,1}.fitter.objective(silja_dat.mrs{1,1},silja_dat.mrs{1,1}.data);
axis([10 3e7 0 max(silja_dat.mrs{1,1}.data.cdf(end,:))]);
[~,~,sleg_h,sleg_text] = legend('Location','BestOutside','FontSize',fontsizes(4));
sleg_sel = [1:3 5:6 8 9:10];
sdata = silja_dat.mrs{1,1}.data;
newsleg_text = cell(1,4*length(sdata.plot_indices));
for i = 1:length(sdata.plot_indices)
    newsleg_text{4*i-3} = ['Data, t = ' num2str(1+sdata.t(sdata.plot_indices(i)))];
    newsleg_text{4*i-2} = ['Model, t = ' num2str(1+sdata.t(sdata.plot_indices(i)))];
    newsleg_text{4*i-1} = ['Inf. Prim, t = ' num2str(1+sdata.t(sdata.plot_indices(i)))];
    newsleg_text{4*i-0} = ['Inf. Meta, t = ' num2str(1+sdata.t(sdata.plot_indices(i)))];
end
legend(sleg_h(sleg_sel),newsleg_text(sleg_sel),'Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('');
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%ssilja_r_cdf',savedir),'-djpeg');

silja_gomp.mrs{1,1}.fitter.objective(silja_gomp.mrs{1,1},silja_gomp.mrs{1,1}.data);
axis([10 3e7 0 max(silja_gomp.mrs{1,1}.data.cdf(end,:))]);
legend('Location','BestOutside','FontSize',fontsizes(4));
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Tumor Size','FontSize',fontsizes(2));
ylabel('Cumulative Number of Tumors','FontSize',fontsizes(3));
title('');
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
print(sprintf('%ssilja_r_gomp_cdf',savedir),'-djpeg');

plotCDF2_vP1(silja_dat.mrs{1,1},silja_gomp.mrs{1,1},silja_dat.mrs{1,1}.data,...
    1:3,[10 3e7 0 max(silja_dat.mrs{1,1}.data.cdf(end,:))],hpic,wpic,...
    'Tumor Size','Cumulative Number of Tumors',fontsizes,legpos4,'auto','auto');

plotMetaOriginTimes_v1_0(silja_dat.mrs{1,1});
cfig = gcf;
cfig.Position = [769.8000 41.8000 766.4000 740.8000];
xlabel('Days Post Inoculation','FontSize',2*fontsizes(2));
ylabel('Number of Metastasis','FontSize',2*fontsizes(3));
title('');
legend('Location','NorthWest','FontSize',2*fontsizes(4));
legend boxoff;
set(gca,'Box','off');
set(gca,'FontSize',2*fontsizes(1));
print(sprintf('%ssilja_r_met',savedir),'-djpeg');

%% 2) plot of immune-suppressed fish data vs immune-competent fish data
%  plot median relative growth male and female together
%  plot male and female separately
%  plot multiple percentiles
plotTumorPercentiles_v1_9a
plotTumorPercentiles_v1_9b_only5e6Inoc
plotTumorPercentiles_v1_9c_2area
plotTumorPercentiles_v1_9d_2area
plotTumorPercentiles_v1_9e_plmm


%% 3) time to recurrence plot
% needs fontsizes and other variables to be in workspace
loadLongRunData_vP1;
print(sprintf('%siwata_recurrence_after_surgery',savedir),'-djpeg');

loadLongRunAnalysis_vP1;
print(sprintf('%siwata_recurrence_time_estimates',savedir),'-djpeg');

%% 4) sample permutation test
% loadPermutationParameters_v3_2; (this generates gender_diff_results.mat) 
gdres = load('mm_form1X_perm_fits/gender_diff_results_save010621.mat');
figure;
histogram(gdres.mfd,25);
hold on;
plot([1 1]*gdres.mfd0,[0 10],'LineWidth',2);
hold off;
xlabel('Gender Difference Measure');
ylabel('Number of Samples');
annotation('textbox', [0.8, 0.695, 0.105, 0.1], 'String', 'Original Sample');
print(sprintf('%ssamp_permution_kr',savedir),'-djpeg');

% calculate p-value
gender_diff_pval = 1-(1+erf(gdres.zdif0/sqrt(2)))/2;
fprintf('Gender Difference p-value: %5.4f\n',gender_diff_pval);
%% 5) inoculation dosage vs size
plotInitialSize_v1_0
print(sprintf('%sdosage_v_init',savedir),'-djpeg');

endt = clock();
duramins = etime(endt,startt)/60;

fprintf('%s took %3.2f minutes\n',mfilename,duramins);

