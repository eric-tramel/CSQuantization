% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% Set dependencies
clear
csq_deps('common','biht','wavelet','srm');

%% Load in data
X = csq_load_data('image','lena.jpg');
imsize = size(X);
N = imsize(1)*imsize(2);

%% Pre-processing
% Normalization
X = X./norm(X(:));

% Dynamic Range
Xrange = [min(X(:)); max(X(:))];

% Rasterization
x = X(:);

% Mean subtraction
Xmu = mean(x);
x = x - Xmu;

%% Experiment Parameters
% Set up wavelet transform
params.L = log2(min(imsize))-1; % Maximum decomposition
params.imsize = imsize;
params.N = N;
% Additional parameters for bivariate shrinkage
params.end_level = params.L - 1;
params.windowsize = 3;
params.lambda = 50;
params.k = round(0.1*params.N);
% Projection parameters
params.block_based = 1;
params.block_dim = [32 32];
params.Nb = (imsize(1)./params.block_dim(1))*(imsize(2)./params.block_dim(2));
params.subrate = 4;
M = round(params.subrate*N);
params.M = M;
% SRM specific parameters
params.blksize = 32;
params.trans_mode = 'BWHT';
% Recovery parameters
params.htol = 0;
params.maxIter = 3000;
params.verbose = 1;
% Side parameters
blockN = params.block_dim(1)*params.block_dim(2);
params.smoothing = @(z) csq_vectorize( im2col(wiener2(col2im(reshape(z,[blockN params.Nb]),params.block_dim,imsize,'distinct'),[3 3]),params.block_dim,'distinct') );


%% Sparse Transform Function Handles
[psi invpsi] = csq_generate_xform('dwt2d',params);
bs_threshold = csq_generate_threshold('bivariate-shrinkage',params);
top_threshold = csq_generate_threshold('top',params);


%% Random Projection Function Handles
[Phi Phi_t] = csq_generate_projection('gaussian',params);

%% Unification
[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

%% Acquisition
y = sign(Phi(x));

%% Recovery
% Recovery parameters
params.ATrans = AT;
% params.threshold = @(z) top_threshold(bs_threshold(z));
params.threshold = bs_threshold;
params.invpsi = invpsi;

% BIHT Recovery
tic
xhat = biht_1d(y,A,params);
biht_time = toc;

%% Evaluation
xhat = xhat + Xmu;  % Return the mean to the result
% mse = sum(norm(xhat-x).^2)./params.N;
d_snr = SNR(xhat,X(:));

csq_printf('Recovered x in %0.2f sec. SNR = %f dB\n',biht_time,d_snr);
csq_printf('Bitrate = %f bpp.\n',M/params.N);

%% Put up image
figure(1);
imagesc(reshape(xhat,imsize)); 
axis image;
colormap(gray);

