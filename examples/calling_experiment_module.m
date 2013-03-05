% calling_experiment_module.m
% 
% Demonstrates how to call the experiment modules.
clear
csq_deps('srm','wavelet','common','experiments','biht','inpaint');

X = csq_load_data('image','lena.jpg');

% Set BIHT-2D experimental parameters
params.block_based = 1;                         % Block acquisition?
params.block_dim = [32 32];                     % Block acq. dimensions
params.htol = 0;                                % Maximum hamming error
params.maxIter = 200;                          % Recovery iterations
params.threshold = 'bivariate-shrinkage';       % Set threshold type
params.lambda = 20;                             % Required B-S parameter
params.xform = 'dwt2d';                         % Sparse Transform
params.projection = 'gaussian';                  % Projection type
params.blksize = 32;                            % Req. SRM parameter
params.trans_mode = 'BWHT';                     % Req. SRM parameter
params.verbose = 1;
params.k = round(0.05*512*512);
params.L = 4;

target_bitrate = 2;

% Call the module
[XF results] = experiment_module_biht2d(X,target_bitrate,params);