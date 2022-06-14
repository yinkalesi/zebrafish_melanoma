%% optimizePermutation_Stage1_v1_0.m
%  Version 4.0
%  Author: Adeyinka Lesi
%  Contributors: Prof David Rumschitzki
%  Date: 8/31/20
% take data from several summary files and create one distribution
%% Version History
%  1.0: from optimizeWithGetSelectedOptima_v4_2

prefix = 'permutations_s2_083120/';
postfix = '_area_vs_count_over_actprim.csv';

% kg_pre kg_exp kr_pre kr_exp ks_pre km_pre ks_exp cc imm_start imm_ramp
ref_parsf = [3.5546,0.79139,0,0,0.16877,0.00016877,0.47077,Inf,0,0];
ref_parsm = [11.828,0.69859,0,0,0.00036929,3.6929e-07,0.89784,Inf,0,0];

% candidate reduction params
kre = [];
f_kge = [ref_parsf(2) kre];
m_kge = [ref_parsm(2) kre];
f_kgp = ref_parsf(1)*exp(-13*(f_kge-ref_parsf(2)));
m_kgp = ref_parsm(1)*exp(-13*(m_kge-ref_parsm(2)));

fac1 = 10;
fac2 = 10;
fac3 = 2;
fac4 = 3;
dif1 = 0.1;
dif2 = 0.3;

conv = getConverter();
stage1_time = 4;
stage2_time = 1;
stage3_time = 1;
stage4_time = 2;

use_para = true;

perm_list = 1:100;
parfor n = 1:length(perm_list)
    
    j = 1;
    start_parsf = ref_parsf;
    start_parsf(1:2) = [f_kgp(j) f_kge(j)];
    start_parsm = ref_parsm;
    start_parsm(1:2) = [m_kgp(j) m_kge(j)];
    
    lowf =  [start_parsf(1)/fac1 0.25 start_parsf(3)/fac1 0.00 start_parsf(5:6)/fac2 0.0 1e7 3.0 .00];
    highf = [start_parsf(1)*fac1 1.0 start_parsf(3)*fac1 1.0  start_parsf(5:6)*fac2  1.0 1e10 6 15];
    lowm =  [start_parsm(1)/fac1 0.25 start_parsm(3)/fac1 0.00 start_parsm(5:6)/fac2 0.0 1e7 3.0 .00];
    highm = [start_parsm(1)*fac1 1.0 start_parsm(3)*fac1 1.0  start_parsm(5:6)*fac2  1.0 1e10 6 15];
    
    lows = [lowf; lowm];
    highs = [highf; highm];
    
    init12_f = [start_parsf(1:2) 0 0 start_parsf(5:8) 0 0];
    init12_m = [start_parsm(1:2) 0 0 start_parsm(5:8) 0 0];
    init12s = {init12_f,init12_m};
    init34_f = start_parsf([3:4 9:10]);
    init34_m = start_parsm([3:4 9:10]);
    init34s = {init34_f,init34_m};
    gender = {'f','m'};
    label = {sprintf('p%0.3i',perm_list(n)),sprintf('p%0.3i',perm_list(n))};
    for ip = 1:2
        gid = gender{ip};
        lab = label{ip};
        input_file = [prefix lab '_' gid '_dist_RAD_zf2017' postfix];
        input_name = [lab '_' gid '_dist_RAD_zf2017'];
        output_name = ['optima_' gid '_' lab '_090420.txt'];
        dotloc = find(output_name=='.',1);
        record_name_rad1 = [output_name(1:dotloc-1) '_kg1_output.txt'];
        record_name_rad2 = [output_name(1:dotloc-1) '_ks1_output.txt'];
        record_name_rad3 = [output_name(1:dotloc-1) '_kg2_output.txt'];
        record_name_rad4 = [output_name(1:dotloc-1) '_ks2_output.txt'];
        
        disp(output_name);
        
        low1 = lows(ip,:);
        high1 = highs(ip,:);
        
        start_pars1 = init12s{ip};
        disp(printParameters_v1_0(start_pars1));
        [pars_rad1,out_rad1] = getSelectedOptima_v2_1(...
            input_file,input_name,start_pars1,[1:2 6:7],low1,high1,...
            conv,1:7,false,stage1_time,0e2,1,use_para);
        dlmwrite(output_name,pars_rad1,'-append');
        % dlmwrite(record_name_rad1,[out_rad1.opt_out.inputs out_rad1.opt_out.outputs'],'-append');
        
        start_pars2 = pars_rad1;
        
        arg2 = struct('ORIGIN_TIME_CALCULATOR',@getOriginSizes_noShed_v1_0);
        [pars_rad2,out_rad2] = getSelectedOptima_v2_1(...
            input_file,input_name,start_pars2,[6:7],low1,high1,...
            conv,1:7,false,stage2_time,0e2,2,use_para,arg2);
        dlmwrite(output_name,pars_rad2,'-append');
        % dlmwrite(record_name_rad2,[out_rad2.opt_out.inputs out_rad2.opt_out.outputs'],'-append');
        
        % define new constraints
        low2 = low1;
        high2 = high1;
        low2(1) = start_pars2(1)/fac3;
        high2(1) = start_pars2(1)*fac3;
        low2(2) = max(start_pars2(2)-dif1,0);
        high2(2) = min(start_pars2(3)+dif1,1);
        low2(5:6) = start_pars2(5:6)/fac4;
        high2(5:6) = start_pars2(5:6)*fac4;
        low2(7) = max(start_pars2(7)-dif2,0);
        high2(7) = min(start_pars2(7)+dif2,1);
        
        start_pars3 = pars_rad2;
        [pars_rad3,out_rad3] = getSelectedOptima_v2_1(...
            input_file,input_name,start_pars3,[1:2],low2,high2,...
            conv,1:7,false,stage3_time,0e2,1,use_para);
        dlmwrite(output_name,pars_rad3,'-append');
        % dlmwrite(record_name_rad3,[out_rad3.opt_out.inputs out_rad3.opt_out.outputs'],'-append');

        start_pars4 = pars_rad3;
        [pars_rad4,out_rad4] = getSelectedOptima_v2_1(...
            input_file,input_name,start_pars4,[6:7],low2,high2,...
            conv,1:7,false,stage4_time,0e2,2,use_para,arg2);
        dlmwrite(output_name,pars_rad4,'-append');
        % dlmwrite(record_name_rad4,[out_rad4.opt_out.inputs out_rad4.opt_out.outputs'],'-append');
    end
end