% 
% function reconstructed_image = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols)
% 
%   This function performs SPL reconstruction of y using a DCT
%   sparsity basis. Phi gives the projection matrix. The reconstructed
%   image, of size num_rows x num_cols, is returned as
%   reconstructed_image.
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


function reconstructed_image = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols)

[M N] = size(Phi);
block_size = sqrt(N);

Psi = DCT2D_Matrix(block_size);

lambda = 6;
TOL = 0.0001;
D_prev = 0;

num_factor = 0;
max_iterations = 200;

x = Phi' * y;

for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi, Psi, ...
      block_size, num_rows, num_cols, lambda);

  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    if num_factor == 4;
      break;
    end
    lambda = lambda * 0.6;
    num_factor = num_factor + 1;
  end
  D_prev = D;
end

reconstructed_image = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distict');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, Psi, ...
    block_size, num_rows, num_cols, lambda)

x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct'); 
x_hat = wiener2(x, [3, 3]);
x_hat = im2col(x_hat, [block_size block_size], 'distinct');

x_hat = x_hat + Phi' * (y - Phi * x_hat);

x1 = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

x_check = Psi' * x_hat;
threshold = lambda * sqrt(2 * log(num_rows * num_cols)) * ...
    (median(abs(x_check(:))) / 0.6745);
x_check(abs(x_check) < threshold) = 0;

x_bar = Psi * x_check;
x = x_bar + Phi' * (y - Phi * x_bar);

x2 = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

D = RMS(x1, x2);

