function key = changeKey_v3_1(defaultKey,changes)
%% changeKey_v3_1
%  Version 3.1
%  Author: Adeyinka Lesi
%  Date: 5/12/20
%  Project: Tumor Growth, Logarithmic Continuum Form
%
%  key = changeKey_v3_1(defaultKey,changes)
%  defaultKey: struct, variables needed for run
%  changes: struct, variables in key needed to change
%% Version History
%  2.0: recalculates time_zero_sizes
%  2.1: add carrying capacity
%  2.2: time_zero_sizes change and death ramp implementation
%  2.3: time_zero_dist provided by time_zero_calculator
%  2.4: change to key.PARAMETERS
%  3.0: updating so it works with CTC version of solver; also using 
%  key.getParameterStruct funciton
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
% key field
key.FIELD = key.FIELD_GENERATOR(key);

% key file prefix
key.FILE_PREFIX = [key.SET_NAME '\' key.FIELD.name '\' key.FIELD.name];



