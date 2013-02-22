
% Set dependencies
clear
csq_deps('common','biht','wavelet','srm','ssim','qamp');

% Load in data
X = csq_load_data('image','lena.jpg');
imsize = size(X);
N = prod(imsize);

% Preprocessing
X = X./norm(X(:));

% Dynamic range
X_range = [min(X(:)); max(X(:));];

% Rasterization
x = X(:);

% Mean subtraction
% Xmu = mean(x);
% x = x - Xmu;

% Experiment parameters
% Wavelet transform
params.L = floor(log2(min(imsize)))-1;
params.imsize = imsize;
params.N = N;

% % additional parameters for bivariate shringkage
% params.end_level = params.L -1;
% params.windowsize = 3;
% params.lambda = 50;
% params.k = round(0.1*params.N);

% Projection parameter
params.block_based = 0;
params.block_dim = [32 32];
params.Nb = round(prod(imsize)/prod(params.block_dim));
params.subrate = 0.7;
M = round(params.subrate*N);
params.M = M;
% SRM specific parameters
params.blksize = 32;
params.trans_mode = 'BWHT';
% Recovery parameters
params.htol = 0.0001;
params.maxIter = 3000;
params.verbose = 1;
% Side parameters
blockN = prod(params.block_dim);
params.smoothing = @(z) ...
  csq_vectorize( im2col(wiener2(col2im(reshape(z,[blockN params.Nb]), ...
  params.block_dim,imsize,'distinct'),[3 3]),params.block_dim,'distinct') );

% Sparse transform functional handle
[psi invpsi] = csq_generate_xform('dwt2d',params);
soft_threshold = csq_generate_threshold('soft',params);

%% Random Projection Function Handles
[Phi Phi_t] = csq_generate_projection('srm-fft',params);

%% Unification
[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

%% Acquisition
% y = sign(Phi(x));
y = Phi(x);

%% Recovery
% Recovery parameters
params.ATrans = AT;
params.Phi = Phi;
params.Phi_t = Phi_t;
% params.threshold = @(z) top_threshold(bs_threshold(z));
params.threshold = soft_threshold;
params.invpsi = invpsi;
params.psi = psi;

% BIHT Recovery
tic
xhat = CS_AMP_Decoder(y,A,params);
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
figure(1);
imagesc(xhat); 
axis image;
colormap(gray);

