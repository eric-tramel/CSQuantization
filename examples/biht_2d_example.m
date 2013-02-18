% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

% Set dependencies
clear
csq_deps('common','biht','wavelet','srm');

% Load in data
X = csq_load_data('image','cameraman.jpg');
imsize = size(X);
x = X(:);

% Normalization
% x = x ./ norm(x);

% Set up wavelet transform
params.L = log2(min(imsize))-1; % Maximum decomposition
params.imsize = imsize;
params.N = imsize(1)*imsize(2);
% Additional parameters for bivariate shrinkage
params.end_level = params.L - 1;
params.windowsize = 3;
params.lambda = 5;
params.k = round(0.05*params.N);

[psi invpsi] = csq_generate_xform('dwt2d',params);

bs_threshold = csq_generate_threshold('bivariate-shrinkage',params);

% Experiment variables
subrate = 0.5;
M = round(subrate*params.N);

% Generate projection and transform
rand_vect = randperm(params.N)';
select_vect = randperm(round(params.N/2)-1)+1;
select_vect = select_vect(1:round(M/2))';
Phi = @(z) fft1d_f(z, select_vect, rand_vect);
Phi_t = @(z) fft1d_t(z, params.N, select_vect, rand_vect);


[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

% Measurement
y = sign(Phi(x));

% Recovery parameters
params.htol = 0;
params.maxIter = 10;
params.ATrans = AT;
params.verbose = 1;
params.threshold = bs_threshold;
params.invpsi = invpsi;

% BIHT Recovery
tic
xhat = biht_1d(y,A,params);
biht_time = toc;

% Evaluation
mse = sum(norm(xhat-x).^2)./params.N;

fprintf('Recovered x in %0.2f sec. MSE = %f\n',biht_time,mse);
fprintf('Bitrate = %f bpp.\n',M/params.N);

%% Put up image
figure(1);
imagesc(reshape(xhat,imsize)); 
axis image;
colormap(gray);

