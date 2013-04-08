% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% Set dependencies
clear
csq_deps('common','biht','wavelet','srm','ssim');

%% Load in data
X = csq_load_data('image','lena.jpg');
imsize = size(X);
N = imsize(1)*imsize(2);

%% Pre-processing
X = X./norm(X(:));					% Normalization
Xrange = [min(X(:)); max(X(:))];	% Dynamic Range
x = X(:); 							% Rasterization
Xmu = mean(x);						% Mean subtraction
x = x - Xmu;

%% Experiment Parameters
params.N = N;
params.block_based = 1;
params.block_dim = [32 32];
params.imsize = size(X);
params.verbose = 1;

params.projection.id = 'srm-blk';
params.projection.blksize = 32;
params.projection.trans_mode = 'BWHT';
params.projection.subrate = 2;

params.transform.id = 'dwt2d';
params.transform.L = log2(min(imsize))-1;	% Maximum decomposition

params.threshold.id = 'bivariate-shrinkage';
params.threshold.lambda = 50;
params.threshold.end_level = params.transform.L - 1;

params.smoothing.id = 'wiener';

params.biht.htol = 0;
params.biht.maxIter = 75;



%% Sparse Transform Function Handles
[psi invpsi] = csq_generate_xform(params);
bs_threshold = csq_generate_threshold(params);

%% Random Projection Function Handles
[Phi Phi_t] = csq_generate_projection(params);

%% Unification
[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

%% Smoothing function
smoothing = csq_generate_smoothing(params);

%% Acquisition
y = sign(Phi(x));
M = length(y) ./ length(x);

%% Recovery
% BIHT Recovery
tic
xhat = biht_1d(y,A,AT,psi,invpsi,bs_threshold,smoothing,params);
biht_time = toc;

%% Attempt Rescaling
xhat = xhat + Xmu;  % Return the mean to the result
xhat = (xhat - min(xhat)) ./ (max(xhat) - min(xhat));
xhat = (Xrange(2) - Xrange(1)).*(xhat + Xrange(1));
xhat = xhat ./ norm(xhat);
xhat = reshape(xhat,imsize);

resc_X = 255*(X - Xrange(1))./(Xrange(2) - Xrange(1));
resc_xhat = 255*(xhat - min(xhat(:))) ./ (max(xhat(:)) - min(xhat(:)));

%% Evaluation
d_mse = MSE(X,xhat);
d_snr = SNR(xhat(:),X(:));
d_rms = RMS(X,xhat);
% d_ssim = ssim(X,xhat,[0.01 0.03],fspecial('gaussian', 11, 1.5),1);
d_ssim = ssim(resc_X,resc_xhat);

csq_printf('Recovered x in %0.2f sec.\n',biht_time);
csq_printf('SNR = %f dB\n',d_snr);
csq_printf('MSE = %f, (Log: %f)\n',d_mse,log10(d_mse));
csq_printf('RMS = %f, (Log: %f)\n',d_rms,log10(d_rms));
csq_printf('SSIM = %f\n',d_ssim);
csq_printf('Bitrate = %f bpp.\n',M/params.N);

%% Put up image
figure(1); clf;
imagesc(xhat); 
axis image;
colormap(gray);

