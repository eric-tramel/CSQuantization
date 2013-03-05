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
  
%     x_hat = params.smoothing(x);
% x = col2im(x, params.block_dim , params.imsize, 'distinct'); 
x = reshape(x, params.imsize);
x_hat = wiener2(x, [3, 3]);
% x_hat = vectorize(im2col(x_hat, params.block_dim, 'distinct'));
% x_hat = vectorize(x_hat);
x_hat = x_hat(:);
    

    x_hat = x_hat + Phi_t((y - Phi(x_hat)));
%     x_hat = x_hat + A_t(y - A(x_hat));
    x_check = Psi(x_hat);
    x_check = params.threshold(x_check);
%     x_check = Psi(params.threshold(Psi_t(x_hat)));
    x_bar = Psi_t(x_check);
    x = x_bar + Phi_t ((y - Phi(x_bar)));
%     x = x_check + A_t(y-A(x_check));
    D = RMS(x_hat, x);
    
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


function x_check = SPLBivariateShrinkage(x_check, end_level, lambda)

windowsize  = 3;
windowfilt = ones(1, windowsize)/windowsize;

tmp = x_check{1}{3};
Nsig = median(abs(tmp(:)))/0.6745;

for scale = 1:end_level
  for dir = 1:3
    Y_coefficient = x_check{scale}{dir};
    
    Y_parent = x_check{scale+1}{dir};
    
    Y_parent = expand(Y_parent);
    
    Wsig = conv2(windowfilt, windowfilt, (Y_coefficient).^2, 'same');
    Ssig = sqrt(max(Wsig-Nsig.^2, eps));
    
    T = sqrt(3)*Nsig^2./Ssig;
    
    x_check{scale}{dir} = bishrink(Y_coefficient, ...
	Y_parent, T*lambda);
  end
end
