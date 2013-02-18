function psnr = run_cspl_ddwt_biht

addpath(genpath('c:\Dropbox\work\Common\BCS-SPL'));
addpath(genpath('c:\Dropbox\work\Common\Sequences'));
addpath(genpath('c:\Dropbox\work\Common\WaveletSoftware'));
addpath(genpath('c:\Dropbox\work\Common\SRM'));
addpath(genpath('../toolboxes/BCS-SPL-1.5-1'));
global original_image

original_image = double(imread('goldhill.pgm'));
[num_rows num_cols] = size(original_image);

subrate = 0.9;

N = num_rows*num_cols;
M = round(N*subrate);

rand_vect = randperm(N)';
select_vect = randperm(round(N/2)-1)+1;
select_vect = select_vect(1:round(M/2))';
Phi = @(z) fft1d_f(z, select_vect, rand_vect);
Phi_t = @(z) fft1d_t(z, N, select_vect, rand_vect);

x = original_image(:);
y = Phi(x);
yq = sign(y);

x_hat = CS_PL_DDWT_Decoder_1bit(yq, Phi, Phi_t, num_rows, num_cols,2);

reconstructed_image = x_hat;
psnr = PSNR(reconstructed_image, original_image);


function reconstructed_image = ...
  CS_PL_DDWT_Decoder_1bit(y, Phi, Phi_t, num_rows, num_cols, num_levels, max_iterations)

lambda = 50;

if (nargin < 7)
  max_iterations = 200;
end

% set level to have maximum wavelet expansion
if (nargin < 6)
  if floor(log2(num_rows)) < floor(log2(num_cols))
    num_levels = floor(log2(num_rows)) - 3;
  else
    num_levels = floor(log2(num_cols)) - 3;
  end
end

TOL = 0.001;
D_prev = 0;

% x = Phi_t(y);
x = zeros(num_rows*num_cols,1);

num_factor = 1;
for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi, Phi_t, num_rows, num_cols, ...
      lambda, num_levels);

    if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
        lambda = lambda/3;
        num_factor = num_factor + 1;
    end
    D_prev = D;
  
    if(num_factor == 4)
        break;
    end
end

reconstructed_image = reshape(x,  [num_rows num_cols]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, Phi_t, num_rows, num_cols, ...
  lambda, num_levels)
global original_image

[Faf, Fsf] = AntonB;
[af, sf] = dualfilt1;

x = x + Phi_t(y - sign(Phi(x)));
x_hat = x + Phi_t(y - sign(Phi(x)));

x1 = reshape(x_hat, [num_rows num_cols]);

x_hat = cplxdual2D(x1, num_levels, Faf, af);

end_level = num_levels -1;
x_hat = SPLBivariateShrinkage(x_hat, end_level, lambda);

x2 = icplxdual2D(x_hat, num_levels, Fsf, sf);
x = x2(:);

D = RMS(x1, x2);

D2 = PSNR(x2,original_image);
disp([num2str(D) '  ' num2str(D2)]);

figure(1); imagesc(x2); colormap gray;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x_check = SPLBivariateShrinkage(x_check, end_level, lambda)

windowsize  = 3;
windowfilt = ones(1, windowsize)/windowsize;

tmp = x_check{1}{1}{1}{1};
Nsig = median(abs(tmp(:)))/0.6745;

for scale = 1:end_level
  for dir = 1:2
    for dir1 = 1:3
      Y_coef_real = x_check{scale}{1}{dir}{dir1};
      Y_coef_imag = x_check{scale}{2}{dir}{dir1};
      Y_parent_real = x_check{scale+1}{1}{dir}{dir1};
      Y_parent_imag = x_check{scale+1}{2}{dir}{dir1};
      Y_parent_real  = expand(Y_parent_real);
      Y_parent_imag  = expand(Y_parent_imag);
      
      Wsig = conv2(windowfilt, windowfilt, (Y_coef_real).^2, 'same');
      Ssig = sqrt(max(Wsig-Nsig.^2, eps));
      
      T = sqrt(3)*Nsig^2./Ssig;
      
      Y_coef = Y_coef_real + sqrt(-1)*Y_coef_imag;
      Y_parent = Y_parent_real + sqrt(-1)*Y_parent_imag;
      Y_coef = bishrink(Y_coef, Y_parent, T*lambda);
      
      x_check{scale}{1}{dir}{dir1} = real(Y_coef);
      x_check{scale}{2}{dir}{dir1} = imag(Y_coef);
    end
  end
end



