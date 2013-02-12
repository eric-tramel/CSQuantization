% 
% function reconstructed_image = ...
%     BCS_SPL_DDWT_Decoder(y, Phi, num_rows, num_cols, num_levels)
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
%   See:
%     S. Mun and J. E. Fowler, "Block Compressed Sensing of Images
%     Using Directional Transforms," submitted to the IEEE
%     International Conference on Image Processing, 2009
%
%   Originally written by SungKwang Mun, Mississippi State University
%

%
% BCS-SPL: Block Compressed Sensing - Smooth Projected Landweber
% Copyright (C) 2009  James E. Fowler
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
    CS_PL_DDWT_Decoder(y, Phi, Phi_t, num_rows, num_cols, nonNegative, num_levels)

  
if nargin < 7
  if floor(log2(num_rows)) < floor(log2(num_cols))
    num_levels = floor(log2(num_rows)) - 3;
  else
    num_levels = floor(log2(num_cols)) - 3;
  end
  if nargin < 6
    nonNegative = 1;
  end
end

lambda = 100; % changed from 25 original BCS-SPL
max_iterations = 1000;

TOL = 0.001;
D_prev = 0;
x = Phi_t(y);

x = reshape(x, [num_rows num_cols]);
row_offset = 0;col_offset = 0;
if num_rows < 2^round(log2(num_rows))
  row_offset = (2^round(log2(num_rows)) - num_rows)/2;
end
if num_cols < 2^round(log2(num_cols))
  col_offset = (2^round(log2(num_cols)) - num_cols)/2;
end

num_factor = 1;
for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi,Phi_t, num_rows, num_cols, ...
      lambda, nonNegative, num_levels, row_offset, col_offset);

    if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
        lambda = lambda/3;
        num_factor = num_factor + 1;
    end
    D_prev = D;
  
    if(num_factor == 5)
        break;
    end
end

reconstructed_image = x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, Phi_t, num_rows, num_cols, ...
    lambda, nonNegative, num_levels, row_offset, col_offset)
global original_image

[Faf, Fsf] = AntonB;
[af, sf] = dualfilt1;

end_level = num_levels - 1;

x_check = symextend_rec(x ,row_offset,col_offset);
x_check = cplxdual2D(x_check,num_levels, Faf, af);
x_check = SPLBivariateShrinkage(x_check, end_level, lambda);
 
x_bar = icplxdual2D(x_check, num_levels, Fsf, sf);
x_bar = x_bar(row_offset+1:row_offset+num_rows,col_offset+1:col_offset+num_cols);
x_bar = x_bar(:);
x = x_bar + Phi_t ((y - Phi(x_bar)));
x = reshape(x, [num_rows num_cols]);

if nonNegative
x(x<0) = 0;
end

D = RMS(x_bar, x(:));
% PSNR( reshape(x,[num_rows num_cols]), original_image) 

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


function y = symextend_rec(x, Nrow, Ncol)
y = [fliplr(x(:,1:Ncol)) x x(:,end:-1:end-Ncol+1)];
y = [flipud(y(1:Nrow,:)); y ;y(end:-1:end-Nrow+1,:)];