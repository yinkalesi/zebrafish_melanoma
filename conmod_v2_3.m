function results = conmod_v2_3(defaultKey,changes)
%% conmod_v2_3
%  Version 2.3
%  Author: Adeyinka Lesi
%  Date: 3/29/18
%  Project: Tumor Growth, Logarithmic Continuum Form
%
%  results = conmod_v2_1(defaultKey,changes)
%  defaultKey: struct, variables needed for run
%  changes: struct, variables in key needed to change
%% Version History
%  1.1: check if FIELD_GENERATOR needs to be used
%  2.0: recalculates time_zero_sizes - necessary for use with
%  getKey_backtrack (also, FILED_GENERATOR always called). Also, cdfs are
%  calculated using prePlotCalcs_v1_2
%  2.1: adding carrying capacity
%  2.2: changes to time_zero_calculator and implementing death ramp
%  2.3: changes to time_zero_calculator; outputing initial tumor dist

key = defaultKey;
changed = fieldnames(changes);
if(~isempty(changed))
    for f = 1:length(changed)
        field = changed{f};
        % check field exists
        if(isfield(key,field))
            key.(field) = changes.(field);
        else
            % erroneous option - throw a warning
            warning([mfilename '-> Option "' field '" not valid']);
        end
    end
end

% need to recalculate some key variables:

% key rates and title
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
                                                
% calculate time_zeros_sizes (need new rates)
if(key.USING_DEATH_STARTER && key.DEATH_START_TIME > 0)
    modkey = key;
    modkey.PARAMETERS.death1 = 0;
    modkey.RATES = modkey.RATE_GENERATOR(modkey);
    [key.TIME_ZERO_SIZES,key.PLACEMENT_DIST,key.PLACEMENT_RADIUS,key.TIME_ZERO_DIST] = key.TIME_ZERO_CALCULATOR(modkey);
else
    [key.TIME_ZERO_SIZES,key.PLACEMENT_DIST,key.PLACEMENT_RADIUS,key.TIME_ZERO_DIST] = key.TIME_ZERO_CALCULATOR(key);
end

% reset transform
key.TRANSFORM = key.TRANSFORMER(key);
% key field (need new time_zero_sizes)
key.FIELD = key.FIELD_GENERATOR(key);

% key file prefix (need field structure)
key.FILE_PREFIX = [key.SET_NAME '\' key.FIELD.name '\' key.FIELD.name];

% CALCULATE RESULTS
results = key.SOLUTION_FUNCTION(key);
results.key = key;
results = prePlotCalcs_v1_2(results); % need results.key

