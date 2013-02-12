% 
% function reconstructed_image = ...
%     BCS_SPL_DDWT_Decoder(y, Phi, num_rows, num_cols, num_levels, ...
%     max_iterations)
% 
%   This function performs SPL reconstruction of y using a DDWT
%   sparsity basis. Phi gives the projection matrix. The reconstructed
%   image, of size num_rows x num_cols, is returned as
%   reconstructed_image.
%
%   The DDWT transform used for reconstruction had num_levels
%   decompositions levels and is implemented using the cplxdual2D
%   program from the Wavelet Software package from
%   http://taco.poly.edu/WaveletSoftware/. The SPL reconstruction uses
%   bivariate shrinkage to perform the thresholding step of the
%   Landweber iteration. The Wavelet Software package must be in
%   Matlab's current search path for both the DDWT as well as
%   bivariate shrinkage.
%
%   max_iterations is an optional parameter specifying the number of
%   iterations of SPL reconstruction to run; default = 200.
% 
%   See:
%     S. Mun and J. E. Fowler, "Block Compressed Sensing of Images
%     Using Directional Transforms," submitted to the IEEE
%     International Conference on Image Processing, 2009
%
%   Originally written by SungKwang Mun, Mississippi State University
%

%
% BCS-SPL: Block Compressed Sensing - Smooth Projected Landweber
% Copyright (C) 2009-2012  James E. Fowler
% http://www.ece.mstate.edu/~fowler
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%


function reconstructed_image = ...
    BCS_SPL_DDWT_Decoder(y, Phi, num_rows, num_cols, num_levels, max_iterations)

lambda = 25;

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

[M N] = size(Phi);
block_size = sqrt(N);

TOL = 0.0001;
D_prev = 0;

x = Phi' * y;

for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
      lambda, num_levels);

  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    break;
  end

  D_prev = D;
end

if (max_iterations > 1)
  [x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
      lambda, num_levels, 'last');
end

reconstructed_image = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distict');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
    lambda, num_levels, last)

if (exist('nor_dualtree.mat'))
  load nor_dualtree
else
  normaliz_coefcalc_dual_tree
end

[Faf, Fsf] = AntonB;
[af, sf] = dualfilt1;

x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct'); 
x_hat = wiener2(x, [3, 3]);
x_hat = im2col(x_hat, [block_size block_size], 'distinct');

x_hat = x_hat + Phi' * (y - Phi * x_hat);

x1 = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

x_check = normcoef(cplxdual2D(symextend(x1, 2^(num_levels - 1)), ...
    num_levels, Faf, af), num_levels, nor);

if (nargin == 9)
  end_level = 1;
else
  end_level = num_levels - 1;
end
x_check = SPLBivariateShrinkage(x_check, end_level, lambda);

x_bar = icplxdual2D(unnormcoef(x_check, num_levels, nor), num_levels, Fsf, sf);
Irow = (2^(num_levels - 1) + 1):(2^(num_levels - 1) + num_rows);
Icol = (2^(num_levels - 1) + 1):(2^(num_levels - 1) + num_cols);
x_bar = x_bar(Irow, Icol);
x_bar = im2col(x_bar, [block_size block_size], 'distinct');

x = x_bar + Phi' * (y - Phi * x_bar);

x2 = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');
  
D = RMS(x1, x2);
  

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
