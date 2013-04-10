function [xhat results] = experiment_module_bcsspl(X,params)
% function results = experiment_module_bcsspl
%
% CSQ Experimental Module for BIHT-2D. For more information regarding the
% CSQ Experimental Module format, please refer to the CSQEM documentation.


% This code demonstrates the usage of the bcsspl code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% THIS CODE IS NOT YET IN AN INTEGRATED AND FUNCTIONAL FORM
error('experiment_module_bcsspl:Unimplemented','This code is not yet working!');

%% Set dependencies
csq_deps('common','wavelet','ssim','qbcsspl','bcs-spl');

%% Setup
csq_required_parameters(params,'subrate', ...
                               'block_based',...
                               'projection', ...
                               'transform', ...
                               'experiment',...
                               'smoothing',...
                               'randseed');

% Make sure that at least the id's are set
csq_required_parameters(params.threshold,'id');
csq_required_parameters(params.projection,'id');
csq_required_parameters(params.transform,'id');
csq_required_parameters(params.smoothing,'id');
csq_required_parameters(params.qbcsspl,'bits','tol','maxIter','quant',...
                                       'meanSubtraction')
% Experimental requirements
csq_required_parameters(params.experiment,'target_bitrate');

% Get image size
imsize = size(X);
params.imsize = imsize;
x = X(:);

% Mean subtraction
if params.qbcsspl.meanSubtraction
  Xmu = mean(x);
else
  Xmu = 0;
end
x = x - Xmu;

if params.block_based == 0
  error('experiment_module_bcsspl:NotBB','BCS-SPL can only run in block-based mode.');
end

% Do we really need this to be a parameter?
% params.Nb = round(prod(imsize)/prod(params.block_dim));

% This section needs to be reworked. This section should not be
% taking in a vector of subrates etc. This would be better implemented
% with an external lookup table function.
% params.subrate = params.subrates(round(target_bitrate*10));
% M = round(params.subrate*params.Nb);
[params.projection.subrate params.qbccspl.bits] = subrate_bit_LUT(params.experiment.target_bitrate);
% ---------------------------------------------

% Should not store the image within the parameters structure. Too much data.
% params.original_image = x;

% Use smoothing generation, here
% params.smoothing = @(z) csq_vectorize(wiener2(reshape(z,imsize),[3 3]));
smoothing = csq_generate_smoothing(params);

% This random seed needs to be updated, perhaps a glue function needs to 
% be made to update the random seed settings without overburdening the 
% experiment module code?
if ~params.randseed, randn('seed',params.randseed); end

% Need to change these settings so that they are defaults which are only
% used if they aren't already specified in the parameters structure
params.L = log2(min(imsize)) - 3; % max L - 3 produces the best
params.end_level = params.L - 1;
params.windowsize = 3;
params.lambda = 10; %ddwt and dwt

switch params.transform.id
  case 'dct2d-blk'
    params.threshold.lambda = 0.3; % 0.8 with block size [8 8] hard thresholding gives good results
                        % 0.3 with block size [8 8] soft thresholding gives good results
    params.threshold.id = 'soft';
    threshold = csq_generate_threshold(params);
    threshold_final = csq_generate_threshold(params);
  case 'dwt2d'
    params.threshold.id = 'bivariate-shrinkage';
    threshold = csq_generate_threshold(params);

    tmp_params = params; 
    tmp_params.threshold.id = 'bivariate-shrinkage-final';
    threshold_final = csq_generate_threshold(tmp_params);
  case 'ddwt2d'
    params.threshold.id = 'ddwt-bs';
    threshold = csq_generate_threshold(params);

    tmp_params = params; 
    tmp_params.threshold.id = 'ddwt-bs-final';
    threshold_final = csq_generate_threshold(tmp_params);
   
  otherwise
    error('unknown type');
end

% Modernized with parameter passing only
[Psi PsiT] = csq_generate_xform(params);
[Phi PhiT] = csq_generate_projection(params);

%% Unification of Phi and Psi
[A AT] = csq_unify_projection(Phi,PhiT,Psi,PsiT);

%% Acquisition
y = Phi(x);
params.projection.M = length(y(:));

switch params.qbcsspl.quant
  case 'sq'
    [yq rate] = SQ_Coding(y, params.qbcsspl.bits, imsize(1),imsize(2));
  case 'dpcm'   
    [yq rate] = DPCM_Coding(y, params.qbcsspl.bits, imsize(1),imsize(2));
  otherwise
    yq = y;
    rate = 0;
end
%% Recovery


% bcsspl_decoder needs to be updated so that it does not
% require these parameters and they should be, instead, passed
% as discrete arguments to the decoder function.
% params.AT = AT;
% params.PsiT = PsiT;
% params.Psi = Psi;

tic
  [xhat iterations] = bcsspl_decoder(yq,A,AT,Psi,PsiT,threshold,smoothing,params);
results.run_time = toc;

%% Attempt Rescaling
xhat = xhat + Xmu;  % Return the mean to the result

%% Finishing
% Outputs
results.iterations = iterations;
results.params = params;
results.true_bitrate = rate;
results.target_bitrate = params.experiment.target_bitrate;






