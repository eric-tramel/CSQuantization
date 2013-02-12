% 
% function reconstructed_image = ...
%     BCS_SPL_DWT_Decoder(y, Phi, num_rows, num_cols, num_levels)
% 
%   This function performs SPL reconstruction of y using a DWT
%   sparsity basis. Phi gives the projection matrix. The reconstructed
%   image, of size num_rows x num_cols, is returned as
%   reconstructed_image.
%
%   The DWT transform used for reconstruction had num_levels
%   decompositions levels and is implemented using the dwt2D
%   program from the Wavelet Software package from
%   http://taco.poly.edu/WaveletSoftware/. The SPL reconstruction uses
%   bivariate shrinkage to perform the thresholding step of the
%   Landweber iteration. The Wavelet Software package must be in
%   Matlab's current search path for both the DWT as well as
%   bivariate shrinkage.
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
    BCS_SPL_DWT_Decoder(y, Phi, num_rows, num_cols, num_levels)

lambda = 20;
max_iterations = 200;

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

end_level = 1;
[x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
    lambda, num_levels, 'last');

reconstructed_image = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distict');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
    lambda, num_levels, last)

[af, sf] = farras;
% [af, sf] = AntonB;
% af = af{1};
% sf = sf{1};
L = 12;

x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct'); 
x_hat = wiener2(x, [3, 3]);
x_hat = im2col(x_hat, [block_size block_size], 'distinct');

x_hat = x_hat + Phi' * (y - Phi * x_hat);

x1 = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

x_check = dwt2D(symextend(x1, L * 2^(num_levels - 1)), ...
    num_levels, af);

if (nargin == 9)
  end_level = 1;
else
  end_level = num_levels - 1;
end
x_check = SPLBivariateShrinkage(x_check, end_level, lambda);

x_bar = idwt2D(x_check, num_levels, sf);
Irow = (L * 2^(num_levels - 1) + 1):(L * 2^(num_levels - 1) + num_rows);
Icol = (L * 2^(num_levels - 1) + 1):(L * 2^(num_levels - 1) + num_cols);
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
