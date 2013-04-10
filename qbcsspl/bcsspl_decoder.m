% 
% function reconstructed_image = ...
%     BCS_SPL_DWT_Decoder(y, Phi, num_rows, num_cols, num_levels)



function reconstructed_image = ...
    bcsspl_decoder(y, A, AT, Psi, PsiT, threshold, smoothing,params)

max_iterations = params.maxIter; %200;
TOL = params.tol; %0.0001;

% Get intial prediction of the recovered image via back-projection.
x = AT(y);

D_prev = 0;
for i = 1:max_iterations
  
 [x D] = spl_recovery(y, x, A, AT, Psi, PsiT, params);
    
  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    break;
  end
  if (params.verbose)
    D_psnr = PSNR(params.original_image, PsiT(x));
    csq_printf('%d: RMS: %.2f, PSNR: %.2f (dB) \n', i, D,D_psnr);
  end
  D_prev = D;
end
params.end_level = 1;
x = spl_recovery(y, x, A, AT, Psi, PsiT, params);

x = PsiT(x);
reconstructed_image = reshape(x,params.imsize);

function [x D] = spl_recovery(y, x, A, AT, Psi, PsiT, params)

x_hat = Psi(params.smoothing(PsiT(x)));
x_hat = x_hat + AT(y - A(x_hat));
if params.end_level~=1
  x_check = params.threshold(x_hat);
else % if final, only use baseband and the respective surrounding subbands.
  x_check = params.threshold_final(x_hat);  
end
x = x_check + AT(y-A(x_check));
D = RMS(x_hat, x);
