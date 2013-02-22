function psnr = run_bcsspl_ddwt_biht

addpath(genpath('c:\Dropbox\work\Common\BCS-SPL'));
addpath(genpath('c:\Dropbox\work\Common\Sequences'));
addpath(genpath('c:\Dropbox\work\Common\WaveletSoftware'));
addpath(genpath('../toolboxes/BCS-SPL-1.5-1'));
global original_image

original_image = double(imread('lenna.pgm'));
[num_rows num_cols] = size(original_image);

subrate = 4;
block_size = 32;
N = block_size*block_size;
M = round(N*subrate);
% Phi = orth(randn(N,M))';
Phi = round(rand(M,N));
Phi(Phi<0.5) = -1;

x = im2col(original_image, [block_size block_size],'distinct');

y = Phi*x;
yq = sign(y);

x_hat = BCS_SPL_DDWT_Decoder_1bit(yq, Phi, num_rows, num_cols);

reconstructed_image = x_hat;
psnr = PSNR(reconstructed_image, original_image);


function reconstructed_image = ...
  BCS_SPL_DDWT_Decoder_1bit(y, Phi, num_rows, num_cols, num_levels, max_iterations)

lambda = 1;

if (nargin < 6)
  max_iterations = 200;
end

% set level to have maximum wavelet expansion
if (nargin < 5)
  if floor(log2(num_rows)) < floor(log2(num_cols))
    num_levels = floor(log2(num_rows)) - 3;
  else
    num_levels = floor(log2(num_cols)) - 3;
  end
end

N = size(Phi,2);
block_size = sqrt(N);

TOL = 0.000001;
D_prev = 0;

% x = Phi' * y;
x = zeros(size(Phi,2),size(y,2));

alpha = norm(Phi'*Phi); % maximum singluar value

% alpha = 1./eye(size(Phi,1),size(Phi,2)).*norm(Phi(1,:));
% alpha = ((Phi*Phi')\Phi)';
num_factor = 1;
for i = 1:max_iterations
  [x D] = SPLIteration(y, x, alpha, Phi, block_size, num_rows, num_cols, ...
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

reconstructed_image = col2im(x, [block_size block_size], ...
  [num_rows num_cols], 'distict');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, alpha, Phi, block_size, num_rows, num_cols, ...
  lambda, num_levels)
global original_image

[Faf, Fsf] = AntonB;
[af, sf] = dualfilt1;

%  x = x + Phi' * (y - sign(Phi * x));
%  x_hat = x + Phi' * (y - sign(Phi * x));
x_hat = x + 1./alpha.*Phi' * (y - sign(Phi * x));
% x_hat = x + alpha* (y - sign(Phi * x));

x1 = col2im(x_hat, [block_size block_size], ...
  [num_rows num_cols], 'distinct');

x_hat = cplxdual2D(x1, num_levels, Faf, af);

end_level = num_levels -1;
x_hat = SPLBivariateShrinkage(x_hat, end_level, lambda);

x2 = icplxdual2D(x_hat, num_levels, Fsf, sf);
x = im2col(x2, [block_size block_size], 'distinct');



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



