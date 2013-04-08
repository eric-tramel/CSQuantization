% This code demonstrates the usage of the bcsspl code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% Set dependencies
clear
csq_deps('common','wavelet','ssim','bcs-spl-dpcm','bcsspl');

%% Load in data
X = csq_load_data('image','lena.jpg');
% X = X(129:384,129:384);
imsize = size(X);
x = X(:);

% type = 'dwt2d';
% type = 'ddwt2d';
type = 'dct2d-blk';

% quant = 'sq';
quant = 'dpcm';
% quant = 'noquant';

% Mean subtraction
% Xmu = mean(x);
Xmu = 0;
x = x - Xmu;

%% Experiment Parameters
params.imsize = imsize;

% Projection parameters
params.block_based = 1;
params.block_dim = [8 8];
params.N = round(prod(imsize)/prod(params.block_dim));

params.subrate = 0.25;
M = round(params.subrate*params.N);
params.M = M;

% Recovery parameters
params.tol = 0.0001;
params.maxIter = 400;
params.verbose = 1;
params.randseed = 0;
params.original_image = x;

params.smoothing = @(z) csq_vectorize(wiener2(reshape(z,imsize),[3 3]));


% paramteres for wavelets
params.L = log2(min(imsize)) - 3; % max L - 3 produces the best
params.end_level = params.L - 1;
params.windowsize = 3;
params.lambda = 10; %ddwt and dwt

switch type
  case 'dct2d-blk'
    params.lambda = 0.3; % 0.8 with block size [8 8] hard thresholding gives good results
                        % 0.3 with block size [8 8] soft thresholding gives good results
    params.threshold = csq_generate_threshold('soft',params);
    params.threshold_final = csq_generate_threshold('soft',params);
  case 'dwt2d'
      params.threshold = csq_generate_threshold('bivariate-shrinkage',params);
      params.threshold_final = csq_generate_threshold('bivariate-shrinkage-final',params);
  case 'ddwt2d'
    params.threshold = csq_generate_threshold('ddwt-bs',params);
    params.threshold_final = csq_generate_threshold('ddwt-bs-final',params);
   
  otherwise
    error('unknown type');
end

[Psi PsiT] = csq_generate_xform(type,params);

% randn('seed',0);
[Phi PhiT] = csq_generate_projection('gaussian',params);

%% Unification of Phi and Psi
[A AT] = csq_unify_projection(Phi,PhiT,Psi,PsiT);

%% Acquisition
y = Phi(x);
B = 5;

switch quant 
  case 'sq'
    [yq rate] = SQ_Coding(y, B, imsize(1),imsize(2));
  case 'dpcm'   
    [yq rate] = DPCM_Coding(y, B, imsize(1),imsize(2));
  otherwise
    yq = y;
    rate = 0;
end
%% Recovery
% Recovery parameters
params.AT = AT;
params.PsiT = PsiT;
params.Psi = Psi;

tic
xhat = bcsspl_decoder(yq,A,params);
recon_time = toc;

%% Attempt Rescaling
xhat = xhat + Xmu;  % Return the mean to the result

%% Evaluation
d_mse = MSE(X,xhat);
d_snr = SNR(xhat(:),X(:));
d_rms = RMS(X,xhat);
% d_ssim = ssim(X,xhat,[0.01 0.03],fspecial('gaussian', 11, 1.5),1);
d_ssim = ssim(X,xhat);
d_psnr = PSNR(xhat,X);

csq_printf('Recovered x in %0.2f sec.\n',recon_time);
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

