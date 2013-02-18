% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

% Set dependencies
clear
csq_deps('common','biht','wavelet','srm');

global invpsi

% Determine sparse transform
xform = 'dwt2d';

X = csq_load_data('image','cameraman.jpg');
x = X(:);
% x = X(:)/norm(X(:));
imsize = size(X);

% Normalization
% x = x ./ norm(x);
N = length(x);

% Set up wavelet transform
xparams.L = 8; %log2(min(size(X))-1
xparams.imsize = imsize;
% Additional parameters for bivariate shrinkage
xparams.end_level = xparams.L - 1;
xparams.windowsize = 3;
xparams.lambda = 5;
xparams.k = round(0.05*N);

[psi invpsi] = csq_generate_xform(xform,xparams);

bs_threshold = csq_generate_threshold('bivariate-shrinkage',xparams);

% Experiment variables
subrate = 0.5;
N = imsize(1)*imsize(2);
M = round(subrate*N);
K = round(0.05*N);

% Generate projection and transform
rand_vect = randperm(N)';
select_vect = randperm(round(N/2)-1)+1;
select_vect = select_vect(1:round(M/2))';
Phi = @(z) fft1d_f(z, select_vect, rand_vect);
Phi_t = @(z) fft1d_t(z, N, select_vect, rand_vect);


[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

% Measurement
y = sign(Phi(x));

% Recovery parameters
params.k = K;
params.htol = 0;
params.maxIter = 400;
params.ATrans = AT;
params.N = N;
params.verbose = 1;
params.threshold = bs_threshold;
params.imsize = size(X);
params.L = xparams.L;

% BIHT Recovery
tic
xhat = biht_1d(y,A,params);
biht_time = toc;

% Evaluation
mse = sum(norm(xhat-x).^2)./N;

fprintf('Recovered x in %0.2f sec. MSE = %f\n',biht_time,mse);
fprintf('Bitrate = %f bpp.\n',M/N);

%% Put up image
figure(1);
imagesc(reshape(invpsi(xhat),imsize)); 
axis image;
colormap(gray);

