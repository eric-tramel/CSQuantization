% 
% function reconstructed_image = ...
%     BCS_SPL_DWT_Decoder(y, Phi, num_rows, num_cols, num_levels)



function reconstructed_image = ...
    bcsspl_decoder(y, A, params)

% lambda = params.lambda; %20;
max_iterations = params.maxIter; %200;
% num_levels = params.L;

% block_size = sqrt(params.Nb);
Phi = params.Phi;
Phi_t = params.Phi_t;
Psi = params.psi;
Psi_t = params.invpsi;

A_t = params.ATrans;

TOL = params.tol; %0.0001;
x = Phi_t(y);

% [af, sf] = farras;

% params.end_levels = params.num_levels-1;
D_prev = 0;
for i = 1:max_iterations
  
    x_hat = params.smoothing(x);
    x_hat = x_hat + Phi_t((y - Phi(x_hat)));
%     x_hat = x_hat + A_t(y - A(x_hat));
    x_check = Psi(x_hat);
    x_check = params.threshold(x_check);
    x_bar = Psi_t(x_check);
    x = x_bar + Phi_t ((y - Phi(x_bar)));
%     x = x_bar + A_t(y-A(x_bar));
    D = RMS(x_hat, x)
    
  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    break;
  end

  D_prev = D;
end

params.end_levels = 1;
x_hat = params.smoothing(x);
x_hat = x_hat + Phi_t((y - Phi(x_hat)));
x_check = Psi(x_hat);
x_check = params.threshold(x_check);
x_bar = Psi_t(x_check);
x = x_bar + Phi_t ((y - Phi(x_bar)));

reconstructed_image = reshape(x,params.imsize);
