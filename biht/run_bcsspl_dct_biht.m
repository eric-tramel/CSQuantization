function psnr = run_bcsspl_dct_biht

% Should use the new csq_deps
%addpath(genpath('c:\Dropbox\work\Common\BCS-SPL'));
%addpath(genpath('../toolboxes/BCS-SPL-1.5-1'));
csq_deps('bcs-spl','common-csq','common-image')

% Using the csq_ glue for loading data
%original_image = double(imread('lenna.pgm'));
original_image = csq_load_data('image','lena.jpg'); % My file location is different !

[num_rows num_cols] = size(original_image);

subrate = 8.0;
block_size = 16;
N = block_size*block_size;
M = round(N*subrate);
% Phi = orth(randn(N,M))';
Phi = round(rand(M,N));
Phi(Phi<0.5) = -1;

x = im2col(original_image, [block_size block_size],'distinct')/norm(original_image(:));

bit = 2;
rate = bit*subrate;
disp(rate);

y = Phi*x;
yq = sign(y);


x_hat = BCS_SPL_DCT_Decoder_1bit(yq, Phi, num_rows, num_cols,original_image);

reconstructed_image = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distict');
psnr = PSNR(reconstructed_image, original_image);

function reconstructed_image = BCS_SPL_DCT_Decoder_1bit(y, Phi, ...
    num_rows, num_cols, original_image)

N = size(Phi, 2);
block_size = sqrt(N);

Psi = DCT2D_Matrix(block_size);

lambda = 0.001; % changing this paramester would increase PSNR
TOL = 0.001;  % changing this paramester would increase PSNR
D_prev = 0;

num_factor = 0;
max_iterations = 2000;

x = zeros(size(Phi,2),size(y,2));
% x = Phi'*y;

for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi, Psi, ...
      block_size, num_rows, num_cols, lambda,original_image);

  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
  if num_factor == 4;
      break;
  end
    lambda = lambda * 0.6;
    num_factor = num_factor + 1;
  end
  D_prev = D;
end

reconstructed_image = x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, Psi, ...
    block_size, num_rows, num_cols, lambda,original_image)

x = x + Phi'*(y - sign(Phi*x));

x_hat = Psi * x;
threshold = lambda * sqrt(2 * log(num_rows*num_cols)) * ...
    (median(abs(x_hat(:))) / 0.6745);
x_hat(abs(x_hat) < threshold) = 0;
x = Psi' * x_hat;

x_bar = x;

x_bar = col2im(x_bar, [block_size block_size], [num_rows num_cols], 'distinct');
x_bar = wiener2(x_bar, [3, 3]);

x = im2col(x_bar, [block_size block_size], 'distinct');
% x = x/norm(x(:)); %normailzation does not work
D = RMS(x_hat, x);
D2 = PSNR(x_bar,original_image);
disp([num2str(D) '  ' num2str(D2)]);

figure(1); imagesc(x_bar); colormap gray;



  
