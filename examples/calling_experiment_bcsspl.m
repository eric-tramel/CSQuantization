% calling_experiment.m
% 
% Demonstrates how to call the experiment modules.
clear

csq_deps('common','experiments','wavelet','ssim','bcs-spl-dpcm','bcsspl');

% Experiment settings
filename = 'lena';
method = 'bcsspl';
xform = 'dwt2d';
quant = 'dpcm';
projection = 'binary';

image = [filename '.jpg'];
bitrates = 0.1:0.1:1.5;

repo_dir = csq_get_repo_dir();
results_filename = [filename '_' method '_' xform '_' quant '_' projection '.mat'];
results_path = [repo_dir '/experiments/results/' results_filename];

% Set bcsspl parameters
params.block_based = 1;          % Block acquisition?
params.tol = 0.0001;
params.maxIter = 400;
params.block_dim = [8 8];        % Block acq. dimensions
params.xform = xform;          % Sparse Transform, 'dwt2d','ddwt2d','dct2d-blk'
params.quant = quant;           % Quanzation method 'sq', 'dpcm', 'noquant'
params.projection = projection;    % Projection type, 'gaussian', 'binary';
params.meanSubtraction = 0;
params.subrates = [0.05 0.08 0.11 0.11 0.14 0.17 0.2 0.25 0.33 0.37 0.4 0.44 0.47 0.51 0.54];
params.bits = [4 5 5 6 6 6 6 6 6 6 6 6 6 6 6 6]; 
params.randseed = 0;                           % fixed random seed, if 0 then random
params.verbose = 1;

% Call the module
image_ratedistortion_experiment(image,bitrates,results_path,@experiment_module_bcsspl,params);
