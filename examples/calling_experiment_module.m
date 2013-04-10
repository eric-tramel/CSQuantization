% calling_experiment_module.m
% 
% Demonstrates how to call the experiment modules.
clear
csq_deps('srm','wavelet','common','experiments','biht','inpaint','proj','bcs-spl');

filename = 'lena.jpg';
X = csq_load_data('image',filename);


% General Parameters
params.rand_seed = 1;							% Seed for the RNG's
params.block_based = 1;                         % Block acquisition?
params.block_dim = [32 32];                     % Block acq. dimensions
params.imsize = size(X);
params.N = length(X(:));
params.verbose = 1;

% CS Projection Parameters
params.projection.id = 'gaussian';
params.projection.blksize = 32;
params.projection.trans_mode = 'BWHT';

% Transform Parameters
params.transform.id = 'dct2d-blk';
params.transform.L = 4;

% Thresholding parameters
params.threshold.id = 'top';
params.threshold.lambda = 30;
params.threshold.k = round(0.05*params.N);
params.threshold.end_level = params.transform.L - 1;

% Smoothing parameters
params.smoothing.id = 'wiener';
params.smoothing.radius = 2;
params.smoothing.window_dim = [3 3];

% BIHT Parameters
params.biht.htol = 0;                                % Maximum hamming error
params.biht.maxIter = 100;                          % Recovery iterations

% Experiment parameters
params.experiment.target_bitrate = 2;

% Call the module
[XF results] = experiment_module_biht2d(X,params);

% Display results
figure(1);
    imagesc(XF);
    axis image;
    colormap(gray);
    xlabel(sprintf('PSNR: %0.2fdB',PSNR(X,XF)));