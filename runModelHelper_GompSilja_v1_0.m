%% runModelHelper_GompSilja_v1_0
%  Version 1.0
%  Author: Adeyinka Lesi
%  Date: 2/17/21
%  Project: Tumor Growth, Logarithmic Continuum Form
%% Version History
%  1.0: looking at Silja's data distributions (made 031218)
%  1.1: looking at later distributions (made 022520) which are divided
%  based on Batch
%  1.2: test of fits

startt = clock();

conv_int = getLinearConverter(1000);
conv_area = getLinearConverter(100);
conv_vol = getConverter();
select1 = 1:3;

save_name = 'Silja_022520_r2_s1_031020_ap_GompSil_021921';

% best so far 0.064169,4.7997e+07,0,0,0.00066628,6.6628e-07,0.89506,Inf,0,0
psil = [0.077512,2.4434e+07,0,0,0.00018252,1.8252e-07,0.98549,Inf,0,0
0.077264,3.2915e+07,0,0,0.00019612,1.9612e-07,0.97136,Inf,0,0
0.060464,3.2915e+07,0,0,0.0010089,1.0089e-06,0.86656,Inf,0,0
0.068595,4.7997e+07,0,0,0.00097533,9.7533e-07,0.87583,Inf,0,0
0.064169,4.7997e+07,0,0,0.00066628,6.6628e-07,0.89506,Inf,0,0
0.083806,2.4913e+07,0,0,0.00062822,6.2822e-07,0.87786,Inf,0,0
0.066616,2.4913e+07,0,0,0.0051647,5.1647e-06,0.74922,Inf,0,0
0.069437,3.9299e+07,0,0,0.0055825,5.5825e-06,0.75375,Inf,0,0
0.068861,3.9299e+07,0,0,0.0052301,5.2301e-06,0.7535,Inf,0,0];

% psil = [2.6158,0.84469,0,0,0.10474,0.00010474,0.50368,Inf,0,0
% 3.0495,0.82416,0,0,0.3669,0.0003669,0.41473,Inf,0,0
% 2.8796,0.83375,0,0,0.4184,0.0004184,0.42188,Inf,0,0
% 3.3372,0.81202,0,0,0.74102,0.00074102,0.36056,Inf,0,0
% 3.4054,0.82193,0,0,0.69585,0.00069585,0.37004,Inf,0,0
% 3.1155,0.81901,0,0,0.75512,0.00075512,0.36308,Inf,0,0
% 3.0877,0.82997,0,0,1.0273,0.0010273,0.34605,Inf,0,0
% 3.0845,0.82464,0,0,0.59085,0.00059085,0.37571,Inf,0,0
% 3.7434,0.80998,0,0,0.66082,0.00066082,0.38209,Inf,0,0
% 3.4423,0.80367,0,0,0.58812,0.00058812,0.38236,Inf,0,0
% 2.7915,0.83545,0,0,0.92107,0.00092107,0.37008,Inf,0,0
% 6.52,0.76783,0,0,0.36631,0.00036631,0.37998,Inf,0,0
% 5.7332,0.75019,0,0,0.037855,3.7855e-05,0.58456,Inf,0,0
% 5.8778,0.77454,0,0,0.027847,2.7847e-05,0.59545,Inf,0,0
% 5.3724,0.75425,0,0,0.79075,0.00079075,0.34788,Inf,0,0
% 2.3883,0.84976,0,0,0.87746,0.00087746,0.36959,Inf,0,0
% 2.4823,0.84639,0,0,0.31199,0.00031199,0.4371,Inf,0,0];


% 
% psil = [7.3816,0.75198,0,0,0.00011519,1.1519e-07,0.97387,Inf,0,0
% 8.359,0.74732,0,0,0.00012319,1.2319e-07,0.97677,Inf,0,0
% 7.3597,0.75222,0,0,0.00014917,1.4917e-07,0.97066,Inf,0,0
% 8.1885,0.7428,0,0,0.0012584,1.2584e-06,0.82159,Inf,0,0
% ];

% psil = [21.437,0.63748,0,0,0.0015172,1.5172e-06,0.80331,Inf,0,0
% 15.611,0.65152,0,0,0.0016284,1.6284e-06,0.81916,Inf,0,0
% 20.728,0.61168,0,0,0.0018675,1.8675e-06,0.81319,Inf,0,0
% 21.066,0.61296,0,0,0.00184,1.84e-06,0.81412,Inf,0,0
% 19.312,0.60556,0,0,0.0014845,1.4845e-06,0.83912,Inf,0,0
% 24.098,0.5953,0,0,0.0055819,5.5819e-06,0.72749,Inf,0,0
% 19.864,0.60953,0,0,0.0051822,5.1822e-06,0.73962,Inf,0,0
% 23.026,0.60186,0,0,0.0064757,6.4757e-06,0.7155,Inf,0,0
% ];

% psil = [12.384,0.71546,0,0,8.5845e-05,8.5845e-08,0.98061,Inf,0,0
% 8.6685,0.74287,0,0,9.7991e-05,9.7991e-08,0.99145,Inf,0,0
% 4.5472,0.8,0,0,7.7352e-05,7.7352e-08,0.99523,Inf,0,0
% 4.4573,0.76562,0,0,0.00014187,1.4187e-07,1,Inf,0,0
% 4.5615,0.8,0,0,6.7169e-05,6.7169e-08,1,Inf,0,0
% 8.0211,0.76525,0,0,0.00014085,1.4085e-07,0.93736,Inf,0,0
% 5.1406,0.78919,0,0,0.00025978,2.5978e-07,0.91253,Inf,0,0
% 4.246,0.77528,0,0,0.00087836,8.7836e-07,0.85667,Inf,0,0
% 4.564,0.8,0,0,7.6288e-05,7.6288e-08,0.9798,Inf,0,0
% 13.608,0.64999,0,0,0.0054859,5.4859e-06,0.7335,Inf,0,0];


datafile1 = 'Silja_022520_r2_s1_031020/dist_RAD_Silja_all_mf_area_vs_count_over_actprim.csv';
datafile2 = 'Silja_022520_r2_s1_031020/dist_RAD_Silja_B1_mf_area_vs_count_over_actprim.csv';
datafile3 = 'Silja_022520_r2_s1_031020/dist_RAD_Silja_B2_mf_area_vs_count_over_actprim.csv';
datafile4 = 'distribution_files_Silja_031218_r3_byT1size/dist_Silja_031218_r3_allT1_area_vs_count_over_prim.csv';
datafile5 = 'distribution_files_Silja_031218_r3_byT1size/dist_Silja_031218_r3_largeT1_area_vs_count_over_prim.csv';
datafiles = {datafile1,datafile2,datafile3,datafile4,datafile5};


mrs = cell(length(datafiles),size(psil,1));
run_list = [1];
if(~exist('objs_func_vals','var'))
    objs_func_vals = zeros(length(datafiles),size(psil,1));
else
    if(size(psil,1)>size(objs_func_vals,2))
        old_objs = objs_func_vals;
        objs_func_vals = zeros(length(datafiles),size(psil,1));
        retain_list = setdiff(1:length(datafiles),run_list);
        objs_func_vals(retain_list,1:size(old_objs,2)) = old_objs(retain_list,:);
    else
        objs_func_vals(run_list,1:size(psil,1)) = 0;
    end
end

Nx = 500;
dt = 0.1;

stage_sel = 1;
tzc_depth = struct('TIME_ZERO_CALCULATOR_DEPTH',2,'ORIGIN_TIME_CALCULATOR',@getOriginSizes_v2_0,'RATE_GENERATOR',@getRates_Gompertz_v2_0);

for rn = run_list
    for i = 1:size(psil,1)
        switch rn
            case 1
                parn = psil(i,:);
                parn([3:4 9:10]) = 0;
                [mrs{rn,i},datan,mg_fitn] = runModelFunction_v5_1(datafiles{rn},parn,conv_vol,Nx,dt,select1,tzc_depth,stage_sel);
                mrs{rn,i}.data = datan;
                mrs{rn,i}.fit = mg_fitn;
                axis([10 1e8 0 max(datan.cdf(end,:))]);
                objs_func_vals(rn,i) = timeWeightedObjective_v1_0(mg_fitn,datan.weights);
%                 plotMetaOriginTimes_v1_0(mrs{rn,i});
%                 plotMetaOriginTimes_v2_0(mrs{rn,i});
                plotMetaOriginTimes_v2_1(mrs{rn,i});
            case 2
                parn = psil(i,:);
                parn([3:4 9:10]) = 0;
                [mrs{rn,i},datan,mg_fitn] = runModelFunction_v5_1(datafiles{rn},parn,conv_vol,Nx,dt,select1,tzc_depth,stage_sel);
                mrs{rn,i}.data = datan;
                mrs{rn,i}.fit = mg_fitn;
                axis([10 1e8 0 max(datan.cdf(end,:))]);
                objs_func_vals(rn,i) = timeWeightedObjective_v1_0(mg_fitn,datan.weights);
                plotMetaOriginTimes_v2_0(mrs{rn,i});
            case 3
                parn = psil(i,:);
                parn([3:4 9:10]) = 0;
                [mrs{rn,i},datan,mg_fitn] = runModelFunction_v5_1(datafiles{rn},parn,conv_vol,Nx,dt,select1,tzc_depth,stage_sel);
                mrs{rn,i}.data = datan;
                mrs{rn,i}.fit = mg_fitn;
                axis([10 1e8 0 max(datan.cdf(end,:))]);
                objs_func_vals(rn,i) = timeWeightedObjective_v1_0(mg_fitn,datan.weights);
                plotMetaOriginTimes_v2_0(mrs{rn,i});
            case 4
                parn = psil(i,:);
                parn([3:4 9:10]) = 0;
                [mrs{rn,i},datan,mg_fitn] = runModelFunction_v5_1(datafiles{rn},parn,conv_vol,Nx,dt,select1,tzc_depth,stage_sel);
                mrs{rn,i}.data = datan;
                mrs{rn,i}.fit = mg_fitn;
                axis([10 1e8 0 max(datan.cdf(end,:))]);
                objs_func_vals(rn,i) = timeWeightedObjective_v1_0(mg_fitn,datan.weights);
                plotMetaOriginTimes_v2_0(mrs{rn,i});
            case 5
                parn = psil(i,:);
                parn([3:4 9:10]) = 0;
                [mrs{rn,i},datan,mg_fitn] = runModelFunction_v5_1(datafiles{rn},parn,conv_vol,Nx,dt,select1,tzc_depth,stage_sel);
                mrs{rn,i}.data = datan;
                mrs{rn,i}.fit = mg_fitn;
                axis([10 1e8 0 max(datan.cdf(end,:))]);
                objs_func_vals(rn,i) = timeWeightedObjective_v1_0(mg_fitn,datan.weights);
                plotMetaOriginTimes_v2_0(mrs{rn,i});
        end
    end
end
save([save_name '_' mfilename],'mrs');

endt = clock();
duramins = etime(endt,startt)/60;

fprintf('%s took %3.2f minutes\n',mfilename,duramins);