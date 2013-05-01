% calling_experiment.m
% 
% Demonstrates how to call the experiment modules.
clear

csq_deps('srm','wavelet','common','experiments','biht','ssim','bcs-spl','proj');

% Experiment settings
image = 'lena.jpg';
bitrates = linspace(0.05,0.99,5);
repo_dir = csq_get_repo_dir();
results_filename = 'lena_biht2d_sub1bpp.mat';
results_path = [repo_dir '/results/' results_filename];

% General Settings
params.rand_seed = 1;
params.block_based = 0;                         % Block acquisition?
params.block_dim = [32 32];                     % Block acq. dimensions
params.verbose = 1;

% Set BIHT-2D parameters
params.biht.htol = 2;                                % Maximum hamming error
params.biht.maxIter = 1000;                          % Recovery iterations

% Projection Parameters
params.projection.id = 'srm-blk';                  % Projection type
params.projection.blksize = 32;                            % Req. SRM parameter
params.projection.trans_mode = 'BWHT';                     % Req. SRM parameter

% Transform Parameters
params.transform.id = 'dwt2d';                         % Sparse Transform
params.transform.L = 4;

% Threshold Parameters
params.threshold.id = 'bivariate-shrinkage';       					% Set threshold type
params.threshold.lambda = 20;                             % Required B-S parameter
params.threshold.k = round(0.05*512*512);

% Smoothing parameters
params.smoothing.id = 'none';
params.smoothing.radius = 2;
params.smoothing.window_dim = [3 3];

% Call the module
image_ratedistortion_experiment(image,bitrates,results_path,@experiment_module_biht2d,params);
