% 
% function reconstructed_image = ...
%     BCS_SPL_DWT_Decoder(y, Phi, num_rows, num_cols, num_levels)



function [reconstructed_image i] = ...
    bcsspl_decoder(y, A, AT, Psi, PsiT, threshold, threshold_final, smoothing, params)

csq_required_parameters(params,'qbcsspl','threshold','imsize');
csq_required_parameters(params.qbcsspl,'maxIter','tol');

% Check to see if we are in some kind of oracle mode
ORACLE = 0;
if isfield(params,'original_image')
  ORACLE = 1;
end

% Check to see if the verbose mode is set
if ~isfield(params,'verbose')
  params.verbose = 1;
end

max_iterations = params.qbcsspl.maxIter; %200;
TOL = params.qbcsspl.tol; %0.0001;

% Get intial prediction of the recovered image via back-projection.
x = AT(y);

D_prev = 0;
for i = 1:max_iterations
  
 [x D] = spl_recovery(y, x, A, AT, Psi, PsiT, threshold, threshold_final, smoothing, params);
    
  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    break;
  end
  if (params.verbose)
    if ORACLE
      D_psnr = PSNR(params.original_image(:), PsiT(x));
      csq_printf('%d: RMS: %.2f, PSNR: %.2f (dB) \n', i, D,D_psnr);
    else
      csq_printf('%d: RMS: %.2f\n', i, D);
    end
  end
  D_prev = D;
end

params.threshold.end_level = 1; % Change the thresholding end level
x = spl_recovery(y, x, A, AT, Psi, PsiT, threshold, threhsold_final, smoothing, params);

x = PsiT(x);
reconstructed_image = reshape(x,params.imsize);

%---------------------------------------------------------------------
function [x D] = spl_recovery(y, x, A, AT, Psi, PsiT, threshold, threshold_final, smoothing, params)

x_hat = Psi(smoothing(PsiT(x)));
x_hat = x_hat + AT(y - A(x_hat));
if params.threshold.end_level ~= 1
  x_check = threshold(x_hat);
else % if final, only use baseband and the respective surrounding subbands.
  x_check = threshold_final(x_hat);  
end
x = x_check + AT(y-A(x_check));
D = RMS(x_hat, x);
