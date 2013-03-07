function [xhat results] = experiment_module_bcsspl(X,target_bitrate,params)
% function results = experiment_module_bcsspl
%
% CSQ Experimental Module for BIHT-2D. For more information regarding the
% CSQ Experimental Module format, please refer to the CSQEM documentation.


% This code demonstrates the usage of the bcsspl code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% Set dependencies
csq_deps('common','wavelet','ssim','bcs-spl-dpcm','bcsspl');

%% Setup
csq_required_parameters(params,'subrates','bits', 'tol', 'maxIter', ...
                           'block_based','projection', 'quant', ...
                           'xform','meanSubtraction', 'randseed');

%% Load in data
% X = csq_load_data('image','lena.jpg');
% X = X(129:384,129:384);

imsize = size(X);
params.imsize = imsize;
x = X(:);

% Mean subtraction
if params.meanSubtraction
  Xmu = mean(x);
else
  Xmu = 0;
end
x = x - Xmu;

%% Experiment Parameters

% Projection parameters
params.block_based = 1;
% params.block_dim = [8 8];
params.Nb = round(prod(imsize)/prod(params.block_dim));

% params.subrate = 0.25;
params.subrate = params.subrates(round(target_bitrate*10));
M = round(params.subrate*params.Nb);
params.M = M;

% Recovery parameters
% params.tol = 0.0001;
% params.maxIter = 400;
% params.verbose = 1;
% params.randseed = 0;
params.original_image = x;

params.smoothing = @(z) csq_vectorize(wiener2(reshape(z,imsize),[3 3]));

if ~params.randseed, randn('seed',params.randseed); end

% paramteres for wavelets
params.L = log2(min(imsize)) - 3; % max L - 3 produces the best
params.end_level = params.L - 1;
params.windowsize = 3;
params.lambda = 10; %ddwt and dwt

switch params.xform
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

[Psi PsiT] = csq_generate_xform(params.xform,params);
[Phi PhiT] = csq_generate_projection(params.projection,params);
%% Unification of Phi and Psi
[A AT] = csq_unify_projection(Phi,PhiT,Psi,PsiT);

%% Acquisition
y = Phi(x);
params.bit_depth = params.bits(round(target_bitrate*10));

switch params.quant
  case 'sq'
    [yq rate] = SQ_Coding(y, params.bit_depth, imsize(1),imsize(2));
  case 'dpcm'   
    [yq rate] = DPCM_Coding(y, params.bit_depth, imsize(1),imsize(2));
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
  [xhat iterations] = bcsspl_decoder(yq,A,params);
results.run_time = toc;

%% Attempt Rescaling
xhat = xhat + Xmu;  % Return the mean to the result

%% Finishing
% Outputs
results.iterations = iterations;
results.params = params;
results.Phi = Phi;
results.Phi_t = PhiT;
results.Psi = Psi;
results.Psi_t = PsiT;
results.true_bitrate = rate;
results.target_bitrate = target_bitrate;






