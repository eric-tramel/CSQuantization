% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

% Set dependencies
clear
csq_deps('common','biht','wavelet','srm');

X = csq_load_data('image','goldhill.jpg');
x = X(:);
imsize = size(X);

% Set up wavelet transform
xparams.L = 5;
xparams.imSize = imsize;

[psi invpsi] = csq_generate_xform('dwt2d',xparams);

% Experiment variables
subrate = 0.5;
N = imsize(1)*imsize(2);
M = round(subrate*N);
K = round(0.25*N);

% Generate projection and transform
rand_vect = randperm(N)';
select_vect = randperm(round(N/2)-1)+1;
select_vect = select_vect(1:round(M/2))';
Phi = @(z) fft1d_f(z, select_vect, rand_vect);
Phi_t = @(z) fft1d_t(z, N, select_vect, rand_vect);


[A AT] = csq_unify_projection(Phi,Phi_t,psi,invpsi);

% Measurement
y = sign(A(x));

% Recovery parameters
params.k = K;
params.htol = 0;
params.maxIter = 3000;
params.ATrans = AT;
params.N = N;
params.verbose = 1;

% BIHT Recovery
tic
xhat = biht_1d(y,A,params);
biht_time = toc;

% Evaluation
mse = sum(norm(xhat-x).^2)./N;

fprintf('Recovered x in %0.2f sec. MSE = %f\n',biht_time,mse);

%% Put up image
figure(1);
imagesc(reshape(xhat,imsize)); colormap(gray);

