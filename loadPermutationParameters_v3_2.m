%% loadPermutationParameters
%  version 3.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 12/8/20
%% Version History
%  load parameters except ones that are exempt
%  2.0: immune-suppressed fish data
%  2.1: immune-competent fish data (zf2017), growth parameters for male and
%  female are the same
%  2.2: immune-competent fish data (zf2017), growth parameters for male and
%  female are different
%  2.3: for paper
%  3.0: for michaelis menten results
%  3.1: form1X results
%  3.2: Michaelis Menten form1X version for paper

startt = clock();

dirname = 'permutations_s2_083120/';
dirname2 = 'ct10_lt015_ht04_021720_rev2_s2_081320/';
fitdirname = 'mm_form1X_perm_fits/';
setname = 'dist_NORAD_zf2017';
prefix = 'optima';
suffix1 = 's3mm1X*_121720.txt';
suffix2 = 's3mm1X*_121720.txt';
descsuff = 'desc.txt';
distsuff = 'area_vs_count_over_actprim.csv';
nperms = 100;
remlist = [];
keeplist = setdiff(1:nperms,remlist);

% load real sample parameters
pfem = [3.3111,0.78887,5.3368,0.00080299,8.8984e-05,8.8984e-08,1,Inf,4.0401,0];
pmal = [15.232,0.67368,14.734,2.7135e-05,7.8204e-05,7.8204e-08,0.99941,Inf,4.2539,0];

% get key with Michaelis-Menten params
fem_form0 = @getRates_michaelisMenten_form123_v1_0;
mal_form0 = @getRates_michaelisMenten_form11_v1_0;
args_fem0 = struct('RATE_GENERATOR',fem_form0);
args_mal0 = struct('RATE_GENERATOR',mal_form0);
conv = getConverter();
fake_data = '0,0.01=1e8;1e8';
% getKeyFromArgs = @(pars,args) runModelFunction_v5_0(fake_data,pars,conv,1000,0.01,[],args,3);
getKeyFromArgs = @(pars,args,datafile) runModelFunction_v5_0(datafile,pars,conv,500,0.01,[1:2 4:6],args,3);

% get kr for original sample
N=100000;
x_range = linspace(1,1e8,N);
datafile_f = sprintf('%s%s_f_%s',dirname2,setname,distsuff);
[fres,datf,fit_f0] = getKeyFromArgs(pfem,args_fem0,datafile_f); close;
obj_f0 = timeWeightedObjective_v1_0(fit_f0,datf.weights);
datafile_m = sprintf('%s%s_m_%s',dirname2,setname,distsuff);
[mres,datm,fit_m0] = getKeyFromArgs(pmal,args_mal0,datafile_m); close;
obj_m0 = timeWeightedObjective_v1_0(fit_m0,datm.weights);
fkr0 = fres.key.RATES.death(x_range);
mkr0 = mres.key.RATES.death(x_range);
min0 = find(fkr0-mkr0<0,1,'last');
min_difx0 = x_range(min0);
% figure;
% loglog(x_range,fkr0,x_range,mkr0,x_range,abs(fkr0-mkr0));

% get male/female difference calculator
mf_dif = @(fkr,mkr) log10(sum(abs(fkr-mkr)));
mfd0 = mf_dif(fkr0,mkr0);

% load permutated parameters
ppm = zeros(nperms,size(pmal,2));
ppf = zeros(nperms,size(pfem,2));

% store male fem diff
mfd = zeros(nperms,1);

% fraction female calculated from desc data
frac_fem = zeros(nperms,1);
frac_mal = zeros(nperms,1);

use_plevel = false;
plevel = 2;
noData = zeros(1,nperms);
min_difx = zeros(1,nperms);
denom_f = zeros(1,nperms);
denom_m = zeros(1,nperms);
% obj_val = zeros(2,nperms);
obj_val = load('mm_form1X_perm_fits/obj_val_010621.csv');
remove_list = []; %[34][11    16    30    34    40    57    83    85]
noData(remove_list) = 1;
perm_list = setdiff(1:nperms,remove_list);
for i = perm_list
    mname = sprintf('%s%s_m_p%03i_%s',fitdirname,prefix,i,suffix2);
    fname = sprintf('%s%s_f_p%03i_%s',fitdirname,prefix,i,suffix1);
    
    mlist = dir(mname);
    mname = [fitdirname mlist(1).name];
    flist = dir(fname);
    fname = [fitdirname flist(1).name];
    
    pim = dlmread(mname);
    pif = dlmread(fname);
    
    if(use_plevel)
        if(size(pim,1) < plevel || size(pif,1) < plevel)
            noData(i) = 1;
        else
            ppm(i,:) = pim(plevel,:);
            ppf(i,:) = pif(plevel,:);
        end
    else
        ppm(i,:) = pim(end,:);
        ppf(i,:) = pif(end,:);
    end
    % load desc file
    mdesc = sprintf('%sp%03i_m_%s_%s',dirname,i,setname,descsuff);
    fdesc = sprintf('%sp%03i_f_%s_%s',dirname,i,setname,descsuff);
    descm = dlmread(mdesc);
    descf = dlmread(fdesc);
    frac_mal(i) = sum(descm(:,3))/size(descm,1);
    frac_fem(i) = sum(descf(:,3))/size(descf,1);
    
    % get powerlaw denominator
    lab = sprintf('p%0.3i',i);
    pw_file_f = ['perm_fits/optima_f_' lab '_s3g2_091520.txt'];
    pw_file_m = ['perm_fits/optima_m_' lab '_s3g2_091520.txt'];
    pw_fit_f = dlmread(pw_file_f);
    pw_fit_m = dlmread(pw_file_m);
    
    denom_f(i) = 1-pw_fit_f(end,4)*0.33/0.45;
    denom_m(i) = 1-pw_fit_m(end,4)*0.33/0.45;
    
    % get formX rates
    fem_formi = @(key) getRates_michaelisMenten_form1X_v1_0(key,denom_f(i));
    mal_formi = @(key) getRates_michaelisMenten_form1X_v1_0(key,denom_m(i));
    args_femi = struct('RATE_GENERATOR',fem_formi);
    args_mali = struct('RATE_GENERATOR',mal_formi);
    
    % calculate male female difference
    %     datafile_fi = sprintf('%sp%03i_f_%s_%s',dirname,i,setname,distsuff);
    %     [fresi,datfi,fit_fi] = getKeyFromArgs(ppf(i,:),args_femi,datafile_fi); close;
    %     obj_val(1,i) = timeWeightedObjective_v1_0(fit_fi,datfi.weights);
    %     datafile_mi = sprintf('%sp%03i_m_%s_%s',dirname,i,setname,distsuff);
    %     [mresi,datmi,fit_mi] = getKeyFromArgs(ppm(i,:),args_mali,datafile_mi); close;
    %     obj_val(2,i) = timeWeightedObjective_v1_0(fit_mi,datmi.weights);
    datafile_fi = fake_data;
    fresi = getKeyFromArgs(ppf(i,:),args_femi,datafile_fi); close;
    datafile_mi = fake_data;
    mresi = getKeyFromArgs(ppm(i,:),args_mali,datafile_mi); close;
    fkri = fresi.key.RATES.death(x_range);
    mkri = mresi.key.RATES.death(x_range);
    mfd(i) = mf_dif(fkri,mkri);
    mini = find(fkri-mkri<0,1,'last');
    if(isempty(mini))
        mini = 1;
    end
    min_difx(i) = x_range(mini);
    
    %     % plot rate curves and difference
    %     if(mod(i,10)==1)
    %         figure('units','normalized','outerposition',[.05 .05 .90 .90]);
    %     end
    %     subplot(2,5,mod(i-1,10)+1);
    %     loglog(x_range,fkri,x_range,mkri,x_range,abs(fkri-mkri));
    %     axis([1 1e8 1 inf]);
    %     title(sprintf('%3.2e',mfd(i)));
    
end

% remove indices with no data
datai = find(noData==0);
ppm = ppm(datai,:);
ppf = ppf(datai,:);
nperms = length(datai);
mfd = mfd(datai);

% calculate normalized differences
mfd_avg = mean(mfd);
mfd_std = std(mfd);
zdif = (mfd-mfd_avg)/mfd_std;
zdif0 = (mfd0-mfd_avg)/mfd_std;

save([fitdirname 'gender_diff_results'],'mfd','mfd_avg','mfd_std','mfd0',...
    'zdif','zdif0','obj_val');

% plot resulting histograms

figure;
histogram(zdif,25)
hold on;
plot([1 1]*zdif0,[0 10],'LineWidth',2);
hold off;
xlabel('Nomalized Gender Difference');
ylabel('Number of Samples');
% annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
print([fitdirname 'gender_diff_norm'],'-djpeg');

figure;
histogram(mfd,25);
hold on;
plot([1 1]*mfd0,[0 10],'LineWidth',2);
hold off;
xlabel('Gender Difference Measure');
ylabel('Number of Samples');
% annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
print([fitdirname 'gender_diff_intd'],'-djpeg');

% figure;
% sel = 3;
% histogram(ppf(:,sel),25);
% hold on;
% plot([1 1]*pfem(sel),[0 10],'LineWidth',2);
% hold off;
% xlabel('Female K1 Value');
% ylabel('Count');
% % annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
% print([fitdirname 'K1_hist_F'],'-djpeg');
% 
% figure;
% sel = 4;
% histogram(ppf(:,sel),25);
% hold on;
% plot([1 1]*pfem(sel),[0 10],'LineWidth',2);
% hold off;
% xlabel('Female K2 Value');
% ylabel('Count');
% % annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
% print([fitdirname 'K2_hist_F'],'-djpeg');

% figure;
% histogram(ppf(:,3)./ppf(:,4),25);
% hold on;
% plot([1 1]*pfem(3)/pfem(4),[0 10],'LineWidth',2);
% hold off;
% xlabel('Female K1/K2');
% ylabel('Count');
% % annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
% print([fitdirname 'K1K2ratio_hist_F'],'-djpeg');


% figure;
% sel = 3;
% histogram(ppm(:,sel),25);
% hold on;
% plot([1 1]*pmal(sel),[0 10],'LineWidth',2);
% hold off;
% xlabel('Male K1 Value');
% ylabel('Count');
% % annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
% print([fitdirname 'K1_hist_M'],'-djpeg');
% 
% figure;
% sel = 4;
% histogram(ppm(:,sel),25);
% hold on;
% plot([1 1]*pmal(sel),[0 10],'LineWidth',2);
% hold off;
% xlabel('Male K2 Value');
% ylabel('Count');
% % annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
% print([fitdirname 'K2_hist_M'],'-djpeg');

% figure;
% histogram(ppm(:,3)./ppm(:,4),25);
% hold on;
% plot([1 1]*pmal(3)/pmal(4),[0 10],'LineWidth',2);
% hold off;
% xlabel('Male K1/K2');
% ylabel('Count');
% % annotation('textbox', [0.70, 0.755, 0.1, 0.1], 'String', 'Original Sample');
% print([fitdirname 'K1K2ratio_hist_M'],'-djpeg');

endt = clock();
duramins = etime(endt,startt)/60;

fprintf('%s took %3.2f minutes\n',mfilename,duramins);




