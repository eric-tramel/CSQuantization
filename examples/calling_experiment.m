% calling_experiment.m
% 
% Demonstrates how to call the experiment modules.
clear

csq_deps('srm','wavelet','common','experiments','biht','ssim');

% Experiment settings
image = 'lena.jpg';
bitrates = linspace(0.75,0.99,10);
repo_dir = csq_get_repo_dir();
results_filename = 'lena_biht2d_sub1bpp.mat';
results_path = [repo_dir '/experiments/results/' results_filename];

% Set BIHT-2D parameters
params.block_based = 0;                         % Block acquisition?
params.block_dim = [32 32];                     % Block acq. dimensions
params.htol = 2;                                % Maximum hamming error
params.maxIter = 4000;                          % Recovery iterations
params.threshold = 'bivariate-shrinkage';       % Set threshold type
params.lambda = 50;                             % Required B-S parameter
params.xform = 'dwt2d';                         % Sparse Transform
params.projection = 'srm-blk';                  % Projection type
params.blksize = 32;                            % Req. SRM parameter
params.trans_mode = 'BWHT';                     % Req. SRM parameter
params.verbose = 1;

% Call the module
image_ratedistortion_experiment(image,bitrates,results_path,@experiment_module_biht2d,params);
