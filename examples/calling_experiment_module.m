% calling_experiment_module.m
% 
% Demonstrates how to call the experiment modules.
clear

csq_deps('srm','wavelet','common','experiments','biht');

X = csq_load_data('image','lena.jpg');

% Set BIHT-2D experimental parameters
params.block_based = 0;                         % Block acquisition?
params.block_dim = [32 32];                     % Block acq. dimensions
params.htol = 0;                                % Maximum hamming error
params.maxIter = 3000;                          % Recovery iterations
params.threshold = 'bivariate-shrinkage';       % Set threshold type
params.lambda = 200;                             % Required B-S parameter
params.xform = 'dwt2d';                         % Sparse Transform
params.projection = 'srm-blk';                  % Projection type
params.blksize = 32;                            % Req. SRM parameter
params.trans_mode = 'BWHT';                     % Req. SRM parameter
params.verbose = 1;

target_bitrate = 0.75;

% Call the module
results = experiment_module_biht2d(X,target_bitrate,params);