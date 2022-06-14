% run optimization on several data files
%% Version History
%  1.5: seeking local optima for parameters, separate starting point for
%  male and female
%  1.6: changed how limits are selected and minor mods
%  2.0: quick mod to optimize newseg data
%  2.1: getSelectedOptima_v1_16
%  2.2: different inital parameters; no carrying capacity
%  2.11: true) - Adjusted MMD
%  2.12: using LLS fitter
%  2.21: fewer grids, smaller dt
%  2.22: try 3 stage fitting
%  2.25: without ramp, log distance weights
%  2.26: without ramp/cc, log distance weights
%  2.27: without ramp/cc, log distance weights for s1 and s3 but no weights
%  for s2
%  2.28: letting slight variation happen in growth for immunity fitting as
%  well
%  3.0: implementing time zero size depth > 1
%  3.3: using metastasis birth fitting
%  3.4: using cdf fit with set meta params to fit CC
%  3.5: reduction parameter fitting
%  4.3: reduction parameters and immunity start time
%  4.4: reduction parameters with start time and ramp

prefix = 'ct10_lt015_ht04_021720_rev2_s2_081320/';
postfix = '_area_vs_count_over_actprim.csv';

% kg_pre kg_exp kr_pre kr_exp ks_pre km_pre ks_exp cc imm_start imm_ramp
ref_parsf = [3.3111,0.78887,100,0.65813,8.8984e-05,8.8984e-08,1,Inf,4.2,3];
ref_parsm = [15.232,0.67368,90,0.65813,7.8204e-05,7.8204e-08,0.99941,Inf,4.2,3];

% candidate reduction params
kre = linspace(0.3,0.6,3);
f_kre = [ref_parsf(4) kre];
m_kre = [ref_parsm(4) kre];
f_krp = ref_parsf(3)*exp(-13*(f_kre-ref_parsf(4)));
m_krp = ref_parsm(3)*exp(-13*(m_kre-ref_parsm(4)));

fac1 = 10;
fac2 = 10;
fac3 = 1.2;
dif1 = 0.02;

conv = getConverter();
stage1_time = 0.5;

use_para = false;

for j = 1:length(f_kre)
    
    start_parsf = ref_parsf;
    start_parsf(3:4) = [f_krp(j) f_kre(j)];
    start_parsm = ref_parsm;
    start_parsm(3:4) = [m_krp(j) m_kre(j)];
    
    lowf =  [start_parsf(1)/fac1 0.25 start_parsf(3)/fac1 0.00 start_parsf(5:6)/fac2 0.0 1e7 3.0 .01];
    highf = [start_parsf(1)*fac1 1.0 start_parsf(3)*fac1 1.0  start_parsf(5:6)*fac2  1.0 1e10 6 7];
    lowm =  [start_parsm(1)/fac1 0.25 start_parsm(3)/fac1 0.00 start_parsm(5:6)/fac2 0.0 1e7 3.0 .01];
    highm = [start_parsm(1)*fac1 1.0 start_parsm(3)*fac1 1.0  start_parsm(5:6)*fac2  1.0 1e10 6 7];
    
    lows = [lowf; lowm];
    highs = [highf; highm];
    
    init12_f = start_parsf;
    init12_m = start_parsm;
    init12s = {init12_f,init12_m};
    init34_f = start_parsf([3:4 9:10]);
    init34_m = start_parsm([3:4 9:10]);
    init34s = {init34_f,init34_m};
    gender = {'f','m'};
    label = {['kre' sprintf('%0.3i',round(start_parsf(4)*100))],['kre' sprintf('%0.3i',round(start_parsm(4)*100))]};
    for ip = 1:2
        gid = gender{ip};
        lab = label{ip};
        output_name = ['optima_' gid '_' lab '_081820_apnw_kr_wocc.txt'];
        dotloc = find(output_name=='.',1);
        record_name_norad = [output_name(1:dotloc-1) '_kr_output.txt'];
        
        disp(output_name);
        
        low1 = lows(ip,:);
        high1 = highs(ip,:);
        
        start_pars1 = init12s{ip};
        disp(printParameters_v1_0(start_pars1));
        
        [pars_norad,out_norad] = getSelectedOptima_v2_3(...
            [prefix 'dist_NORAD_zf2017_' gid postfix],['dist_NORAD_zf2017_' gid],...
            start_pars1,[3:4 9:10],low1,high1,conv,[1:2 4:6],true,stage1_time,0e2,3,use_para);
        dlmwrite(output_name,pars_norad,'-append');
        dlmwrite(record_name_norad,[out_norad.opt_out.inputs out_norad.opt_out.outputs'],'-append');
    end
end