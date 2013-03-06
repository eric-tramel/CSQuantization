function x = CS_AMP_Decoder(y, Ain, params)

% Variable Initialization
% Default Variables
htol = 0;
maxIter = 1000;
verbose = 0;
conv = 1e-10;

% Flags
FLAG_Nspecified = 0;
FLAG_ATransspecified = 0;
FLAG_Kspecified = 0;

% Check the input parameters
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

% if isfield(params,'invpsi')
%     invpsi = params.invpsi;
% else
%     invpsi = @(x) x;
% end


% Input handling
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
    Phi = @(x) params.Phi(x);
    Phi_t = @(x) params.Phi_t(x);
    Psi = @(x)  params.psi(x);
    Psi_t = @(x)  params.invpsi(x);
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

% Recovery
% Initialization
x = zeros(N,1); 
% x = AT(y);
hiter = Inf;
iter = 0;
conv_check = Inf;
M = length(y);

norm_y = norm(y);
z = y;
% Main Recovery Loop
% while (htol < hiter) && (iter < maxIter) && (conv_check > conv)

for i = 1:maxIter
    % Convergence checking
    xprev = x;
    
    x_check = Psi(x);
    gamma = x_check + AT(z); %Psi(Phi_t(z))
    tau = largestElement(abs(gamma(:,1)), M);
    x_check = soft_threshold(gamma, tau);
    
    temp = Psi_t(etaprime(gamma(:,1), tau));
    temp = (z/M)*sum(temp(:));
    z = y - A(x_check) + temp; %Phi*Psi'*x_check
    
    x = Psi_t(x_check);
%     disp(PSNR(original_image, col2im(x, [block_size block_size], ...
%        [num_rows num_cols], 'distinct')));

    % Normalize
%     x = x ./ norm(x);
    
    % Evaluate
    hiter = nnz(y - A(x));
    
    if verbose
        csq_printf('[biht_1d.m] Iter %d: \\delta = %f, hiter = %d, ||g||_2 = %f.\n',iter,nnz(x)./N,hiter,norm(gamma));
    end
    
    iter = iter + 1;

    conv_check = norm(x-xprev)./N;

     
  if norm(y - Phi(x))/norm_y < htol
      break;
  end

%     Debug code for the two-dimensional case
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
 
    figure(2);
    subplot(1,2,1);
    imagesc(reshape(Psi_t(x),params.imsize));
    axis image;
    subplot(1,2,2);
    imagesc(reshape(abs(Psi_t(gamma)),params.imsize));
    axis image;
    colormap(gray);
%     refresh;
end

% Finishing
x = Psi_t(x);
x = x ./ norm(x); 

if verbose
    csq_printf('[biht_1d.m] Compelted Recovery. Iters = %d, hfinal = %d.\n',iter,hiter);
end

