% calling_experiment_module.m
% 
% Demonstrates how to call the experiment modules.
clear
csq_deps('srm','wavelet','common','experiments','biht','inpaint','proj');

X = csq_load_data('image','lena.jpg');

% Set BIHT-2D experimental parameters
params.rand_seed = 1;							% Seed for the RNG's
params.block_based = 1;                         % Block acquisition?
params.block_dim = [32 32];                     % Block acq. dimensions
params.htol = 0;                                % Maximum hamming error
params.maxIter = 200;                          % Recovery iterations
params.threshold = 'bivariate-shrinkage';       % Set threshold type
params.lambda = 30;                             % Required B-S parameter
params.xform = 'dwt2d';                         % Sparse Transform
params.projection = 'gaussian';                  % Projection type
params.blksize = 32;                            % Req. SRM parameter
params.trans_mode = 'BWHT';                     % Req. SRM parameter
params.verbose = 1;
params.k = round(0.05*512*512);
params.L = 4;
params.smooth_id = 'deblock';					% Name of smoothing function
params.radius = 2;								% Radius for deblocking filter
params.window_dim = [3 3];						% Weiner filter window size

target_bitrate = 2;

% Call the module
[XF results] = experiment_module_biht2d(X,target_bitrate,params);

% Display results
figure(1);
    imagesc(XF);
    axis image;
    colormap(gray);
    xlabel(sprintf('PSNR: %0.2fdB',PSNR(X,XF)));