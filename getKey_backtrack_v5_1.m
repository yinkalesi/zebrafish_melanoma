function key = getKey_backtrack_v5_1(pars,tf,N_t,data,N_x,x_ref,lowx,unk_dist,name)
%% getKey_backtrack_v5_1
%  Version 5.1
%  Author: Adeyinka Lesi
%  Date: 5/12/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%  
%  key = getKey_unitspaced_v1_0(kg,kr,ks,km,tf,N,init,name)
%  kg: 1x2 vector; growth parameters
%  kr: 1x2 vector; death parameters
%  ks: scalar; shedding parameter
%  km: 1x2 vector; metastasis parameters
%  tf: scalar; simulation end time
%  N_t: scalar; number of time points
%  x_init: 1xM vector; sizes of tumors initially on grid
%  d_init: 1xM vector; number of tumors at each size
%  N_x: scalar; number of size intervals on grid (number of x points - 1)
%  x_ref: scalar; reference size for primary tumor at t=0
%  lowa: scalar; threshold size for data
%  name: string; set name
%% Version History
%  1.0: allow backtracking to take into account the metastasis generation
%  between start time and first observation
%  1.1: 4/25/17 made part where time_zero_sizes are calculated into separate
%  function (no new version due to laziness)
%  1.1: allow backtracking to take into account the metastasis generation between start time and first observation 
%  2.0: implement carrying capacity
%  3.0: if tumor growth rate doesn't make sense for time zero sizes, will not
%  put that particular tumor in init_dist; implemented ramp function for
%  death rate
%  4.0: changing how time_zero_sizes are calculated (checking t2
%  distribution as well to see if some tumors from there are similarly
%  undetected)
%  4.1: adding death_ramp params to parameters list
%  5.0: Allows use of getTransformedDistribution_CTC_v1_0, which allows
%  explicit CTC calculation and reseeding (requiring three new parameters)
%  5.1: using getTimeZeroSizes_v4_1

key = struct();

% helper functions used:
% generate rates from constants
key.RATE_GENERATOR = @getRates_powerLaw_v5_0; 

% solver used
key.SOLUTION_FUNCTION = @getTransformedDistribution_CTC_v1_0;

% plotting function used
key.PLOTTING_FUNCTION = @plotResults_v1_0;

% saving function
key.SAVE_FUNCTION = @saveResults_v1_0;

% variable transormation functions
key.TRANSFORMER = @getTransformationRule_v1_4;

% field setup function
key.FIELD_GENERATOR = @setUpField_v8_0;

% data convertion function (if raw data is not in tumor cell count)
key.CONV = data.conv;

% get parameters from pars
kg = pars(1:2);
kr = pars(3:4);
ks = pars([5 7]);
if(length(pars) < 8)
    % no carrying capacity, shed_meta_ratio, or reseeding parameters
    cc = inf;
    kc = [pars(6) inf 0 0];
    if(pars(6) ~= 0)
        kc(2) = pars(5)/pars(6);
    end
    ramp_center = 0;
    ramp_width = 0;
elseif(length(pars) < 11)
    % no shed_meta_ratio or reseeding parameters
    cc = pars(8);
    kc = [pars(6) inf 0 0];
    if(pars(6) ~= 0)
        kc(2) = pars(5)/pars(6);
    end
    ramp_center = pars(9);
    ramp_width = pars(10);
elseif(length(pars) < 12)
    % no reseeding parameters
    cc = pars(8);
    kc = [pars(6) pars(11) 0 0];
    ramp_center = pars(9);
    ramp_width = pars(10);
else
    cc = pars(8);
    kc = pars([6 11:13]);
    ramp_center = pars(9);
    ramp_width = pars(10);
end
% growth, death and metastasis parameters
key.GROWTH_PARAMETER1 = kg(1);
key.GROWTH_PARAMETER2 = kg(2);
key.DEATH_PARAMETER1 = kr(1);
key.DEATH_PARAMETER2 = kr(2);
key.SHED_PARAMETER1 = ks(1);
key.SHED_PARAMETER2 = ks(2);
key.CARRYING_CAPACITY = cc;
key.CTC_PARAMETER1 = kc(1);
key.SHED_META_RATIO = kc(2);
key.CTC_PARAMETER2 = (key.SHED_META_RATIO-1)*key.CTC_PARAMETER1;
key.SEED_PARAMETER1 = kc(3);
key.SEED_PARAMETER2 = kc(4);

% time parameters
key.TIME_LIMIT = tf;
key.TIME_STEP = tf/N_t;

% define x_init, d_init
x_init = data.x(data.selector{1});
d_init = data.dist(data.selector{1},1);

% placing tumors onto grid
key.INITIAL_TUMOR_DIST = d_init;
key.INITIAL_TUMOR_SIZES = x_init;
key.INITIAL_NUMBER_OF_TUMORS = sum(key.INITIAL_TUMOR_DIST);
% key.HARD_TUMOR_SIZE_LIMIT = max(x_init)*5+2; % placeholder (replaced below)
% key.NUMBER_OF_SIZE_INTERVALS = N_x;
key.DETECTION_LIMIT = lowx;
% def: placement distribution
unk_dist_key = struct('uniform',0,'Gaussian',1,'wedge',2,'Poisson',3);
key.PLACEMENT_DIST = zeros(size(x_init));
% def: placement radius - half-width or 6*std
key.PLACEMENT_RADIUS = .05*ones(size(x_init));

% naming parameters
key.SET_NAME = name;
key.CASE = 1;

% model buffer parameters
% the model only solves up to the point where there is a measurable
% distribution value
key.SOFT_LIMIT_BUFFER = 100;
key.LIMIT_RESET_THRESHOLD = 1E-80;

% death cutoff parameters
% death cutoff is where the death parameter becomes zero after a given
% cutoff size
key.USING_DEATH_CUTOFF = false;
key.DEATH_RATE_CUTOFF_SIZE = 0;

% surgery parameters
% after the given surgery time, the distribution at a size greater than the
% cutoff is set to zero
key.CONDUCTING_SURGERY = false;
key.SURGERY_TIME = 0;
key.SURGERY_CUTOFF_SIZE = 0;

% death starter parameters
% this allows the model to run with a death parameter of zero until a given
% time, after which the parameter is used
key.USING_DEATH_STARTER = false;
key.DEATH_START_TIME = 0;
% death ramp parameters
key.USING_DEATH_RAMP = ramp_center>0;
key.DEATH_RAMP_CENTER = ramp_center;
key.DEATH_RAMP_WIDTH = ramp_width;
key.DEATH_RAMP_CREATOR = @getLinearRampFunction_v1_0;

% parameters for metastsis generation
key.USING_METAGEN = 0;
key.METAGEN_START_SIZE = 1;
key.METAGEN_END_SIZE = x_ref;
key.METAGEN_MAX_TIME = min(1e4,max(ceil((x_ref^(1-kg(2))-1)/(1-kg(2))/kg(1)),tf));

% parameters for CTC concentration
key.INITIAL_CTC_CONC = 0;
key.USING_EXPLICIT_CTC = length(pars)>10;

% want to have a consistent time
key.TIME_SHIFT = 0;

% time points to save;
key.TIMES_TO_STORE = tf;
key.SAVE_NUMBER = min(round(10^(2+mod(log10(tf),1))),N_t)+1;

% setting for reporting percentage completion
key.PERCENT_COMPLETE_INCREMENT = 0.02;

% get rates and title string
key.getParameterStruct = @(key) struct('growth1',key.GROWTH_PARAMETER1,...
                        'growth2',key.GROWTH_PARAMETER2,...
                        'death1',key.DEATH_PARAMETER1,...
                        'death2', key.DEATH_PARAMETER2,...
                        'shed1',key.SHED_PARAMETER1,...
                        'shed2',key.SHED_PARAMETER2,...
                        'meta1',key.SHED_PARAMETER1/key.SHED_META_RATIO,...
                        'meta2',key.SHED_PARAMETER2,...
                        'shed_meta_ratio',key.SHED_META_RATIO,...
                        'ctc1',key.CTC_PARAMETER1,...
                        'ctc2',key.CTC_PARAMETER2,...
                        'seed1',key.SEED_PARAMETER1,...
                        'seed2',key.SEED_PARAMETER2,...
                        'carrying_capacity',key.CARRYING_CAPACITY,...
                        'death_ramp_center',key.DEATH_RAMP_CENTER,...
                        'death_ramp_width',key.DEATH_RAMP_WIDTH,...
                        'death_cutoff',key.DEATH_RATE_CUTOFF_SIZE,...
                        'death_start_time',key.DEATH_START_TIME);
key.PARAMETERS = key.getParameterStruct(key);
                    
[key.RATES, key.TITLE] = key.RATE_GENERATOR(key);
% get death ramp
key.DEATH_RAMP_FUNCTION = key.DEATH_RAMP_CREATOR(key.DEATH_RAMP_CENTER,...
    key.DEATH_RAMP_WIDTH);

% set real hard tumor size limit
size_vals = [getSizeTrajectory_v3_3(key,max(x_init),tf/1000,tf) data.x']*2+2;
key.HARD_TUMOR_SIZE_LIMIT = ceil(max(size_vals(~isinf(size_vals))));
key.NUMBER_OF_SIZE_INTERVALS = min(N_x,key.HARD_TUMOR_SIZE_LIMIT-1);

% get TIME_ZERO_SIZES using RK4
% need to handle tumors that are known to be present but are below
% detection limit. That implies that the specified size in x_init is
% irrelevant for these 'invisible' tumors and they will be placed based on
% a putative distribution
key.TIME_ZERO_CALCULATOR_DEPTH = 3; % length(data.t);
key.ORIGIN_TIME_CALCULATOR = @getOriginSizes_v1_0;
key.unk_dist_key = unk_dist_key;
key.unk_dist = unk_dist;
key.data = data;
key.TIME_ZERO_UPDATER = @setTimeZeroSizes;
key = key.TIME_ZERO_UPDATER(key);

% transformation of x axis
key.TRANSFORM = key.TRANSFORMER(key);
                                       
% def: field is a struct containing variables describing simulation
% conditions - it uses the values specified above
key.FIELD = key.FIELD_GENERATOR(key);       

key.FILE_PREFIX = [key.SET_NAME '\' key.FIELD.name '\' key.FIELD.name];

key.QUIET_MODE = 0;
end

function key = setTimeZeroSizes(old_key)
key = old_key;
key.TIME_ZERO_CALCULATOR = @(new_key) new_key.ORIGIN_TIME_CALCULATOR(new_key,...
    new_key.data,new_key.TIME_ZERO_CALCULATOR_DEPTH);

if(key.USING_DEATH_STARTER && key.DEATH_START_TIME > 0)
    modkey = key;
    modkey.PARAMETERS.death1 = 0;
    modkey.RATES = modkey.RATE_GENERATOR(modkey.PARAMETERS);
    [key.TIME_ZERO_SIZES,key.PLACEMENT_DIST,key.PLACEMENT_RADIUS,key.TIME_ZERO_DIST,...
        key.TIME_ZERO_ORIGIN,key.META_ORIGIN,key.inferred_data,key.SIZE_PROJECTIONS] = key.TIME_ZERO_CALCULATOR(modkey);
else
    [key.TIME_ZERO_SIZES,key.PLACEMENT_DIST,key.PLACEMENT_RADIUS,key.TIME_ZERO_DIST,...
        key.TIME_ZERO_ORIGIN,key.META_ORIGIN,key.inferred_data,key.SIZE_PROJECTIONS] = key.TIME_ZERO_CALCULATOR(key);
end
end