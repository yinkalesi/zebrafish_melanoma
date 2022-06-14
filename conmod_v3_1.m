function results = conmod_v3_1(defaultKey,changes)
%% conmod_v3_1
%  Version 3.1
%  Author: Adeyinka Lesi
%  Date: 5/12/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%
%  results = conmod_v3_1(defaultKey,changes)
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
%  2.4: re-did parameters
%  3.0: Allows use of getTransformedDistribution_CTC_v1_0, which allows
%  explicit CTC calculation and reseeding (requiring three new parameters)
%  3.1: added key.TIME_ZERO_UPDATER function

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
key.PARAMETERS = key.getParameterStruct(key);
                    
[key.RATES, key.TITLE] = key.RATE_GENERATOR(key);  
% get death ramp
key.DEATH_RAMP_FUNCTION = key.DEATH_RAMP_CREATOR(key.DEATH_RAMP_CENTER,...
                                                    key.DEATH_RAMP_WIDTH);

% calculate time_zeros_sizes
% unk_dist_key = struct('uniform',0,'Gaussian',1,'wedge',2,'Poisson',3);
key = key.TIME_ZERO_UPDATER(key);

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

