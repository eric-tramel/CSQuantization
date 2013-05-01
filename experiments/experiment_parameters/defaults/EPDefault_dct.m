% Default parameter settings for using DCT sparse transform
csq_deps('common');

% Transform Parameters
params.transform.id = 'dct';                         	% Use DCT

% Threshold Parameters
params.threshold.id = 'hard';       					% Set threshold type
params.threshold.lambda = 20;                           % Required BS parameter

