% optimizeWithGetSelectedOptima_mm_v1_0
% run optimization on several data files
%% Version History
%  1.0: from optimizeWithGetSelectedOptima_v4_4 but using Michaelis-Mentin
%  parameters (form11)
%  1.1: Michaelis-Mentin form 113
%  1.1: Michaelis-Mentin form 123

prefix = 'ct10_lt015_ht04_021720_rev2_s2_081320/';
postfix = '_area_vs_count_over_actprim.csv';

% kg_pre kg_exp kr_pre kr_exp ks_pre km_pre ks_exp cc imm_start imm_ramp
ref_parsf = [3.3111,0.78887,4.288,4e-4,8.8984e-05,8.8984e-08,1,Inf,4.2687,0];
ref_parsm = [15.232,0.67368,18.104,4e-3,7.8204e-05,7.8204e-08,0.99941,Inf,4.2908,0];

% candidate reduction params
kre = [];
f_kr1 = [ref_parsf(4) kre];
m_kr2 = [ref_parsm(4) kre];
f_krp = ref_parsf(3)/ref_parsf(4)*f_kr1;
m_krp = ref_parsm(3)/ref_parsm(4)*m_kr2;

fac1 = 10;
fac2 = 10;
fac3 = 5;
dif1 = 0.02;

conv = getConverter();
stage1_time = 5;

use_para = false;
mm_args = struct('RATE_GENERATOR',@getRates_michaelisMenten_form123_v1_0);

for j = 1:length(f_kr1)
    
    start_parsf = ref_parsf;
    start_parsf(3:4) = [f_krp(j) f_kr1(j)];
    start_parsm = ref_parsm;
    start_parsm(3:4) = [m_krp(j) m_kr2(j)];
    
    lowf =  [start_parsf(1)/fac1 0.25 start_parsf(3)/fac1 start_parsf(4)/fac3 start_parsf(5:6)/fac2 0.0 1e7 3.0 .01];
    highf = [start_parsf(1)*fac1 1.0 start_parsf(3)*fac1 start_parsf(4)*fac3  start_parsf(5:6)*fac2  1.0 1e10 6 7];
    lowm =  [start_parsm(1)/fac1 0.25 start_parsm(3)/fac1 start_parsm(4)/fac3 start_parsm(5:6)/fac2 0.0 1e7 3.0 .01];
    highm = [start_parsm(1)*fac1 1.0 start_parsm(3)*fac1 start_parsm(4)*fac3  start_parsm(5:6)*fac2  1.0 1e10 6 7];
    
    lows = [lowf; lowm];
    highs = [highf; highm];
    
    init12_f = start_parsf;
    init12_m = start_parsm;
    init12s = {init12_f,init12_m};
    init34_f = start_parsf([3:4 9:10]);
    init34_m = start_parsm([3:4 9:10]);
    init34s = {init34_f,init34_m};
    gender = {'f','m'};
    label = {['kr1' sprintf('%0.3i',round(start_parsf(3)*100))],['kr1' sprintf('%0.3i',round(start_parsm(3)*100))]};
    for ip = 1:2
        gid = gender{ip};
        lab = label{ip};
        output_name = ['optima_' gid '_' lab '_111720_kr_mm123.txt'];
        dotloc = find(output_name=='.',1);
        record_name_norad = [output_name(1:dotloc-1) '_kr_mm123_output.txt'];
        
        disp(output_name);
        
        low1 = lows(ip,:);
        high1 = highs(ip,:);
        
        start_pars1 = init12s{ip};
        disp(printParameters_v1_0(start_pars1));
        
        [pars_norad,out_norad] = getSelectedOptima_mm_v1_0(...
            [prefix 'dist_NORAD_zf2017_' gid postfix],['dist_NORAD_zf2017_' gid],...
            start_pars1,[3:4 9],low1,high1,conv,[1:2 4:8],true,stage1_time,0e2,3,use_para,mm_args);
        dlmwrite(output_name,pars_norad,'-append');
        dlmwrite(record_name_norad,[out_norad.opt_out.inputs out_norad.opt_out.outputs'],'-append');
    end
end