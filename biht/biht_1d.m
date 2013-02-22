function x = biht_1d(y,Ain,params)
% biht_1d(y,Phi,params)
% Code to implement the basic 1D functionality of the the
% 1-Bit CS recovery procedure: Binary Iterated Hard Thresholding.
% Note: This code is only intended 
%
% Inputs:
%       y		--		A Mx1 vector of 1-bit quantized CS measurements
%       A		--		A MxN projection matrix used to calculate y = Ax,
%                       or a function handle which maps N->M.
%       params	-- 		Structure lgorithm parameters:
%                       params.htol: Integer[0,N] (Default = 0).
%                       params.k: Integer[0,N] (Default = M/4).
%                       params.ATrans: NxM matrix function handle. Optional
%                                      if A is a real matrix.
%                       params.maxIter: Integer[1,...] (Default = 1000)
%                       params.verbose
%
% This code is based on the work presented in 
% L. Jauqes, J. Laska, P. T. Boufouros, and R. G. Baraniuk,
% "Robust 1-Bit Compressive Sensing via Binary Stable Embeddings"
% and in its accompanying demonstration code.


%% Variable Initialization
% Default Variables
htol = 0;
maxIter = 1000;
verbose = 0;
conv = 1e-10;

% Flags
FLAG_Nspecified = 0;
FLAG_ATransspecified = 0;
FLAG_Kspecified = 0;

%% Check the input parameters
if isfield(params,'htol')
	htol = params.htol;
end

if isfield(params,'k')
    FLAG_Kspecified = 1;
end

if isfield(params,'maxIter')
	maxIter = params.maxIter;
end

if isfield(params,'ATrans')
	FLAG_ATransspecified = 1;
end

if isfield(params,'N')
	N = params.N;
    FLAG_Nspecified = 1;
end

if isfield(params,'verbose')
    verbose = params.verbose;
end

if isfield(params,'threshold')
    threshold = params.threshold;
else
    if ~FLAG_Kspecified
        error('biht1d:NoK','Need to know how many coefficients to retain on each iteration.');
    end
    threshold = csq_generate_threshold('top',params);
end

if isfield(params,'invpsi')
    invpsi = params.invpsi;
else
    invpsi = @(x) x;
end


%% Input handling
if isa(Ain,'function_handle')
	A = @(x) sign(Ain(x));
    
    if ~FLAG_Nspecified
        error('biht_1d:MissingParameter','Original dimensionality, N, not specified.');
    end
else
	A = @(x) sign(Ain*x);
    Nn = size(Ain,2);
    if ~FLAG_Nspecified
        N = Nn;
    else
        if N ~= Nn
            error('biht_1d:ParameterMismatch','Projection N and parameter N do not match.');
        end
    end
end

if (FLAG_ATransspecified)
	if isa(params.ATrans,'function_handle')
		AT = @(x) params.ATrans(x);
	else
		AT = @(x) params.ATrans*x;
	end
else
	if isa(Ain,'function_handle')
		error('biht_1d:NoTranspose','No projection transpose function specified.');
	else
		AT = @(x) Ain'*x;
	end
end

%% Recovery
% Initialization
% x = zeros(N,1); 
x = AT(y);
hiter = Inf;
iter = 0;
conv_check = Inf;
M = length(y);

% Main Recovery Loop
while (htol < hiter) && (iter < maxIter) && (conv_check > conv)
    % Convergence checking
    xprev = x;
    
    % Gradient
    g = AT(y - sign(A(x)));
    
    % Step
    r = x + (1/M) .* g;
    
    % Update 
    x = params.smoothing(x);
    x = threshold(r); 
    
    % Normalize
    x = x ./ norm(x);
    
    % Evaluate
    hiter = nnz(y - A(x));
    
    if verbose
        csq_printf('[biht_1d.m] Iter %d: \\delta = %f, hiter = %d, ||g||_2 = %f.\n',iter,nnz(x)./N,hiter,norm(g));
    end
    
    iter = iter + 1;

    conv_check = norm(x-xprev)./N;
    
    %% Debug code for the two-dimensional case
%     figure(1);
%     plot(abs(x)); axis tight; grid on; box on;
%     figure(1);
%     L = params.L;
%     for l=1:L
%         W = csq_dwt_vec2cell(x,params.imsize(1),params.imsize(2),L);
%         subplot(L+1,1,l);
%         imagesc(abs(cell2mat(W{l}))); colormap(jet); axis image;
%     end
%     subplot(L+1,1,L+1);
%     imagesc(abs(W{L+1})); axis image;
%  
%     figure(2);
%     subplot(1,2,1);
%     imagesc(reshape(invpsi(x),params.imsize));
%     axis image;
%     subplot(1,2,2);
%     imagesc(reshape(abs(invpsi(g)),params.imsize));
%     axis image;
%     colormap(gray);
%     refresh;
end

% Finishing
x = invpsi(x);
x = x ./ norm(x); 

if verbose
    csq_printf('[biht_1d.m] Compelted Recovery. Iters = %d, hfinal = %d.\n',iter,hiter);
end












