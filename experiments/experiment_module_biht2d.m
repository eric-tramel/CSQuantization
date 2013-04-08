function [XF results] = experiment_module_biht2d(X,params)
% function results = experiment_module_biht2d
%
% CSQ Experimental Module for BIHT-2D. For more information regarding the
% CSQ Experimental Module format, please refer to the CSQEM documentation.

%% Setup
% Check for basic settings
csq_required_parameters(params,'block_based');
% Check for parameter sets
csq_required_parameters(params,'biht','threshold','projection','transform','smoothing');
% Check for BIHT settings
csq_required_parameters(params.biht,'htol','maxIter');
% Make sure that at least the id's are set
csq_required_parameters(params.threshold,'id');
csq_required_parameters(params.projection,'id');
csq_required_parameters(params.transform,'id');
csq_required_parameters(params.smoothing,'id');
% Experimental requirements
csq_required_parameters(params.experiment,'filename','target_bitrate');

% Assume X is coming in as an image matrix
params.imsize = size(X);
params.N = length(X(:));

% Set wavelet decomposition levels
if ~isfield(params.transform,'L')
    params.transform.L = log2(min(params.imsize))-1;     % Wavelet decomposition levels
end

% Some input checking for different experiment modes
if params.block_based
    csq_required_parameters(params,'block_dim');
end

switch params.threshold
    case 'bivariate-shrinkage'
        csq_required_parameters(params,'lambda');
        params.threshold.end_level = params.L - 1;
        params.threhsold.windowsize = 3;
    case 'top'
        csq_required_parameters(params,'k');
end

switch params.projection
    case 'srm-blk'
        csq_required_parameters(params,'blksize','trans_mode');
end

% Determine subrate from bitrate
%   For the BIHT, since all measurements are 1 bit, the target bitrate (in
%   bpp) uniquely determines the subrate we should use. 
params.projection.subrate = params.experiment.target_bitrate;
params.projection.M = round(params.projection.subrate*params.N); % Get number of measurements

% Generate projection set
[Phi Phi_t] = csq_generate_projection(params);

% Generate transform set
[Psi Psi_t] = csq_generate_xform(params);

% Unification
[A AT] = csq_unify_projection(Phi,Phi_t,Psi,Psi_t);

% Generate threshold
threshold = csq_generate_threshold(params);

% Generate smoothing
smoothing = csq_generate_smoothing(params);

%% Experiment
% Normalization and mean subtraction
Xeng = norm(X(:));
xn = X(:) ./ Xeng;
xmean = mean(xn);
xn = xn - xmean;

% Projection
y = sign(Phi(xn));

% Recovery
tic 
    [XF iterations] = biht_1d(y,A,AT,Psi,Psi_t,threshold,smoothing,params);
results.run_time = toc;

% Adding the mean back
XF = XF + xmean;

% Returning the energy
XF = XF .* Xeng;

% Reshape
XF = reshape(XF,params.imsize);

%% Finishing
% Outputs
results.iterations = iterations;
results.params = params;
% results.Phi = Phi;                % Can be generated from params
% results.Phi_t = Phi_t;            % Can be generated from params
% results.Psi = Psi;                % Can be generated from params
% results.Psi_t = Psi_t;            % Can be generated from params
results.true_bitrate = length(y) ./ params.N;
results.target_bitrate = params.experiment.target_bitrate;



