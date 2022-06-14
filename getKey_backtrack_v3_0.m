function key = getKey_backtrack_v3_0(kg,kr,ks,km,ts,tf,N_t,x_init,d_init,N_x,x_ref,lowx,unk_dist,name)
%% getKey_backtrack_v3_0
%  Version 3.0
%  Author: Adeyinka Lesi
%  Date: 10/09/17
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

key = struct();

% helper functions used:
% generate rates from constants
key.RATE_GENERATOR = @getRates_powerLaw_v3_0; 

% solver used
key.SOLUTION_FUNCTION = @getTransformedDistribution_timeSensitive_v1_0;

% plotting function used
key.PLOTTING_FUNCTION = @plotResults_v1_0;

% saving function
key.SAVE_FUNCTION = @saveResults_v1_0;

% variable transormation functions
key.TRANSFORMER = @getTransformationRule_v1_3;

% field setup function
key.FIELD_GENERATOR = @setUpField_v7_2;

% data convertion function (if raw data is not in tumor cell count)
key.CONV = getNullConverter();

% growth, death and metastasis parameters
key.GROWTH_PARAMETER1 = kg(1);
key.DEATH_PARAMETER1 = kr(1);
key.SHED_PARAMETER1 = ks(1);
if(km(1) ~= 0)
    key.SHED_META_RATIO = ks(1)/km(1);
else
    key.SHED_META_RATIO = 0;
end
key.META_PARAMETER1 = km(1);
key.GROWTH_PARAMETER2 = kg(2);
key.DEATH_PARAMETER2 = kr(2);
key.META_PARAMETER2 = km(2);
key.CARRYING_CAPACITY = inf;

% time parameters
key.TIME_LIMIT = tf;
key.TIME_STEP = tf/N_t;

% placing tumors onto grid
key.INITIAL_TUMOR_DIST = d_init;
key.INITIAL_TUMOR_SIZES = x_init;
key.INITIAL_NUMBER_OF_TUMORS = sum(key.INITIAL_TUMOR_DIST);
key.HARD_TUMOR_SIZE_LIMIT = max(x_init)*5+2;
key.NUMBER_OF_SIZE_INTERVALS = N_x;
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
key.USING_DEATH_RAMP = false;
key.DEATH_RAMP_CENTER = 0;
key.DEATH_RAMP_WIDTH = 0;
key.DEATH_RAMP_CREATOR = @getLinearRampFunction_v1_0;

% parameters for metastsis generation
key.USING_METAGEN = 0;
key.METAGEN_START_SIZE = 1;
key.METAGEN_END_SIZE = x_ref;
key.METAGEN_MAX_TIME = max(ceil((x_ref^(1-kg(2))-1)/(1-kg(2))/kg(1)),tf);

% want to have a consistent time
key.TIME_SHIFT = 0;

% time points to save;
key.TIMES_TO_STORE = tf;
key.SAVE_NUMBER = min(round(10^(2+mod(log10(tf),1))),N_t)+1;

% setting for reporting percentage completion
key.PERCENT_COMPLETE_INCREMENT = 0.02;

% transformation of x axis
key.TRANSFORM = key.TRANSFORMER(key);

% get rates and title string
key.PARAMETERS = struct('growth1',key.GROWTH_PARAMETER1,...
                        'growth2',key.GROWTH_PARAMETER2,...
                        'death1',key.DEATH_PARAMETER1,...
                        'death2', key.DEATH_PARAMETER2,...
                        'shed1',key.SHED_PARAMETER1,...
                        'meta1',key.META_PARAMETER1,...
                        'meta2',key.META_PARAMETER2,...
                        'death_cutoff_present',key.USING_DEATH_CUTOFF,...
                        'death_cutoff',key.DEATH_RATE_CUTOFF_SIZE,...
                        'carrying_capacity',key.CARRYING_CAPACITY);
                    
[key.RATES, key.TITLE] = key.RATE_GENERATOR(key);
% get death ramp
key.DEATH_RAMP_FUNCTION = key.DEATH_RAMP_CREATOR(key.DEATH_RAMP_CENTER,...
                                                    key.DEATH_RAMP_WIDTH);

% get TIME_ZERO_SIZES using RK4
% need to handle tumors that are known to be present but are below
% detection limit. That implies that the specified size in x_init is
% irrelevant for these 'invisible' tumors and they will be placed based on
% a putative distribution
key.TIME_ZERO_CALCULATOR = @(new_key) getTimeZeroSizes_v2_0(new_key,x_init,ts,lowx,unk_dist_key,unk_dist);

if(key.USING_DEATH_STARTER && key.DEATH_START_TIME > 0)
    modkey = key;
    modkey.PARAMETERS.death1 = 0;
    modkey.RATES = modkey.RATE_GENERATOR(modkey.PARAMETERS);
    [key.TIME_ZERO_SIZES,key.PLACEMENT_DIST,key.PLACEMENT_RADIUS,key.DIST_SELECTOR] = key.TIME_ZERO_CALCULATOR(modkey);
else
    [key.TIME_ZERO_SIZES,key.PLACEMENT_DIST,key.PLACEMENT_RADIUS,key.DIST_SELECTOR] = key.TIME_ZERO_CALCULATOR(key);
end
    
% def: field is a struct containing variables describing simulation
% conditions - it uses the values specified above
key.FIELD = key.FIELD_GENERATOR(key);       

key.FILE_PREFIX = [key.SET_NAME '\' key.FIELD.name '\' key.FIELD.name];

key.QUIET_MODE = 0;
