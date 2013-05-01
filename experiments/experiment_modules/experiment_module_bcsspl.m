function [xhat results] = experiment_module_bcsspl(X,params)
% function results = experiment_module_bcsspl
%
% CSQ Experimental Module for BIHT-2D. For more information regarding the
% CSQ Experimental Module format, please refer to the CSQEM documentation.


% This code demonstrates the usage of the bcsspl code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

%% THIS CODE IS NOT YET IN AN INTEGRATED AND FUNCTIONAL FORM
% error('experiment_module_bcsspl:Unimplemented','This code is not yet working!');

%% Set dependencies
csq_deps('common','wavelet','ssim','qbcsspl','bcs-spl');

%% Setup
csq_required_parameters(params,...
                               'block_based',...
                               'projection', ...
                               'transform', ...
                               'experiment',...
                               'smoothing',...
                               'rand_seed');

% Make sure that at least the id's are set
csq_required_parameters(params.threshold,'id');
csq_required_parameters(params.projection,'id');
csq_required_parameters(params.transform,'id');
csq_required_parameters(params.smoothing,'id');
csq_required_parameters(params.qbcsspl,'tol','maxIter','quant',...
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

% Is bit-depth specified?
BITS_SPECIFIED = 0;
if isfield(params.qbcsspl,'bits')
  BITS_SPECIFIED = 1;
end

SUBRATE_SPECIFIED = 0;
if isfield(params.projection,'subrate')
  SUBRATE_SPECIFIED = 1;
end

% If we don't know anything about subrate and bits, look it up from the table to match the bitrate
if ~BITS_SPECIFIED && ~SUBRATE_SPECIFIED
  [params.projection.subrate params.qbcsspl.bits] = subrate_bit_LUT(params.experiment.target_bitrate);
end

if BITS_SPECIFIED && ~SUBRATE_SPECIFIED
  % Figure out what subrate is needed to match the bitrate with a given bitdepth
  params.projection.subrate = params.experiment.target_bitrate / params.qbcsspl.bits;
end

if SUBRATE_SPECIFIED && ~BITS_SPECIFIED
  % Figure out what bitdepth is needed to match the bitrate with a given subrate
  params.qbcsspl.bits = floor(params.experiment.target_bitrate / params.projection.subrate);
  params.qbcsspl.bits = max(params.qbcsspl.bits, 1); % Make sure we don't have a bit-depth of 0
end
% ---------------------------------------------

% Should not store the image within the parameters structure. Too much data.
% params.original_image = x;

% Use smoothing generation, here
% params.smoothing = @(z) csq_vectorize(wiener2(reshape(z,imsize),[3 3]));
smoothing = csq_generate_smoothing(params);

% This random seed needs to be updated, perhaps a glue function needs to 
% be made to update the random seed settings without overburdening the 
% experiment module code?
% if ~params.randseed, randn('seed',params.randseed); end

% Need to change these settings so that they are defaults which are only
% used if they aren't already specified in the parameters structure
params.transform.L = log2(min(imsize)) - 3; % max L - 3 produces the best
params.threshold.end_level = params.transform.L - 1;
params.threshold.windowsize = 3;
params.threshold.lambda = 10; %ddwt and dwt

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
fprintf('Q-BCS-SPL, Bit Depth = %d, Subrate = %0.2f\n',params.qbcsspl.bits,params.projection.subrate);
fprintf('           Target Rate = %0.3f, True Rate = %0.3f\n',params.experiment.target_bitrate,rate);
tic
  [xhat iterations] = bcsspl_decoder(yq,A,AT,Psi,PsiT,threshold,threshold_final,smoothing,params);
results.run_time = toc;

%% Attempt Rescaling
xhat = xhat + Xmu;  % Return the mean to the result

%% Finishing
% Outputs
results.iterations = iterations;
results.params = params;
results.true_bitrate = rate;
results.target_bitrate = params.experiment.target_bitrate;






