function results = experiment_module_biht2d(X,target_bitrate,params)
% function results = experiment_module_biht2d
%
% CSQ Experimental Module for BIHT-2D. For more information regarding the
% CSQ Experimental Module format, please refer to the CSQEM documentation.

%% Setup
csq_required_parameters(params,'htol','maxIter','threshold',...
                           'block_based','projection',...
                           'xform');

% Assume X is coming in as an image matrix
params.imsize = size(X);
params.N = length(X(:));
params.L = log2(min(params.imsize))-1;     % Wavelet decomposition levels
params.smoothing = @(z) z;

% Some input checking for different experiment modes
if params.block_based
    csq_required_parameters(params,'block_dim');
    % Determine the number of blocks in the image
    params.Nb = (imsize(1)./params.block_dim(1))*(imsize(2)./params.block_dim(2));
    % Set smoothing function
    blockN = params.block_dim(1)*params.block_dim(2);
    params.smoothing = @(z) csq_vectorize( im2col(wiener2(col2im(reshape(z,[blockN params.Nb]),params.block_dim,imsize,'distinct'),[3 3]),params.block_dim,'distinct') );
end

switch params.threshold
    case 'bivariate-shrinkage'
        csq_required_parameters(params,'lambda');
        params.end_level = params.L - 1;
        params.windowsize = 3;
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
params.subrate = target_bitrate;
params.M = round(params.subrate*params.N); % Get number of measurements

% Generate projection set
[Phi Phi_t] = csq_generate_projection(params.projection,params);

% Generate transform set
[Psi Psi_t] = csq_generate_xform(params.xform,params);

% Unification
[A AT] = csq_unify_projection(Phi,Phi_t,Psi,Psi_t);

% Generate threshold
params.threshold = csq_generate_threshold(params.threshold,params);

% BIHT Recovery parameters
params.ATrans = AT;
params.invpsi = Psi_t;

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
    XF = biht_1d(y,A,params);
results.run_time = toc;

% Adding the mean back
XF = XF + xmean;

% Returning the energy
XF = XF .* Xeng;

%% Finishing
% Outputs
results.params = params;
results.Phi = Phi;
results.Phi_t = Phi_t;
results.Psi = Psi;
results.Psi_t = Psi_t;
results.XF = reshape(XF,params.imsize);
results.true_bitrate = length(y) ./ params.N;
results.target_bitrate = target_bitrate;



