% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% Set dependencies
clear
csq_deps('common','wavelet','ssim','bcs-spl-dpcm');

%% Load in data
X = csq_load_data('image','lena.jpg');
imsize = size(X);
% N = imsize(1)*imsize(2);

%% Pre-processing
% Normalization
% X = X./norm(X(:));

% Dynamic Range
% Xrange = [min(X(:)); max(X(:))];

% Rasterization
x = X(:);

% Mean subtraction
% Xmu = mean(x);
Xmu = 0;
x = x - Xmu;

%% Experiment Parameters
% Set up wavelet transform
% params.L = log2(min(imsize))-1; % Maximum decomposition
params.L = 5;
params.imsize = imsize;
% params.N = N;
% Additional parameters for bivariate shrinkage
params.end_level = params.L - 1;
params.windowsize = 3;
params.lambda = 20;
% params.k = round(0.1*params.N);
% Projection parameters
params.block_based = 1;
params.block_dim = [32 32];
params.Nb = round(prod(imsize)/prod(params.block_dim));
params.subrate = 0.3;
M = round(params.subrate*params.Nb);
params.M = M;
% Recovery parameters
params.tol = 0.0001;
params.maxIter = 200;
params.verbose = 1;
% Side parameters
blockN = params.block_dim(1)*params.block_dim(2);
params.smoothing = @(z) csq_vectorize( im2col(wiener2(col2im(reshape(z,[blockN params.Nb]),params.block_dim,imsize,'distinct'),[3 3]),params.block_dim,'distinct') );


%% Sparse Transform Function Handles
[psi invpsi] = csq_generate_xform('dwt2d',params);
bs_threshold = csq_generate_threshold('bivariate-shrinkage',params);
% top_threshold = csq_generate_threshold('top',params);

randn('seed',0);
%% Random Projection Function Handles
[Phi Phi_t] = csq_generate_projection('gaussian',params);

%% Unification
[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

%% Acquisition
y = Phi(x);
B = 5;
[yq rate] = DPCM_Coding(y, B, imsize(1),imsize(2));

%% Recovery
% Recovery parameters
params.ATrans = AT;
% params.threshold = @(z) top_threshold(bs_threshold(z));
params.threshold = bs_threshold;
params.invpsi = invpsi;
params.psi = psi;
params.Phi = Phi;
params.Phi_t = Phi_t;

% BIHT Recovery
tic
xhat = bcsspl_decoder(y,A,params);
biht_time = toc;

%% Attempt Rescaling
xhat = xhat + Xmu;  % Return the mean to the result
% xhat = (xhat - min(xhat)) ./ (max(xhat) - min(xhat));
% xhat = (Xrange(2) - Xrange(1)).*(xhat + Xrange(1));
% xhat = xhat ./ norm(xhat);
% xhat = reshape(xhat,imsize);

% resc_X = 255*(X - Xrange(1))./(Xrange(2) - Xrange(1));
% resc_xhat = 255*(xhat - min(xhat(:))) ./ (max(xhat(:)) - min(xhat(:)));

%% Evaluation
d_mse = MSE(X,xhat);
d_snr = SNR(xhat(:),X(:));
d_rms = RMS(X,xhat);
% d_ssim = ssim(X,xhat,[0.01 0.03],fspecial('gaussian', 11, 1.5),1);
d_ssim = ssim(X,xhat);
d_psnr = PSNR(xhat,X);

csq_printf('Recovered x in %0.2f sec.\n',biht_time);
csq_printf('SNR = %f dB\n',d_snr);
csq_printf('PSNR = %f dB\n',d_psnr);
csq_printf('MSE = %f, (Log: %f)\n',d_mse,log10(d_mse));
csq_printf('RMS = %f, (Log: %f)\n',d_rms,log10(d_rms));
csq_printf('SSIM = %f\n',d_ssim);
csq_printf('Bitrate = %f bpp.\n',rate);

%% Put up image
figure(1);
imagesc(xhat); 
axis image;
colormap(gray);

