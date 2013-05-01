% Default parameter settings for using DWT2D sparse transform
csq_deps('common','wavelet');

% Transform Parameters
params.transform.id = 'dwt2d';                         	% Use 2D DWT
params.transform.L = 4;									% Specify decomposition level

% Threshold Parameters
params.threshold.id = 'bivariate-shrinkage';       		% Set threshold type
params.threshold.lambda = 20;                           % Required B-S parameter

