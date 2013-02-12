% 
% function reconstructed_image = ...
%     BCS_SPL_CT_Decoder(y, Phi, num_rows, num_cols, contourlet)
% 
%   This function performs SPL reconstruction of y using a contourlet
%   sparsity basis. Phi gives the projection matrix. The reconstructed
%   image, of size num_rows x num_cols, is returned as
%   reconstructed_image.
%
%   The contourlet transform used for reconstruction is a structure;
%   e.g.:
%     contourlet.nlevels = [4 4 4 4];
%     contourlet.pfilter = '9-7';
%     contourlet.dfilter = 'pkva';
%   These parameters are in the syntax required by the pdfbdec program
%   of the Contourlet Toolbox (http://www.ifp.illinois.edu/~minhdo/software/)
%   The Contourlet Toolbox must appear in Matlab's current search path.
%
%   The SPL reconstruction uses bivariate shrinkage to perform the
%   thresholding step of the Landweber iteration. The Wavelet Software
%   package from http://taco.poly.edu/WaveletSoftware/ must be in
%   Matlab's current search path.

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
    BCS_SPL_CT_Decoder(y, Phi, num_rows, num_cols, contourlet)

lambda = 10;
max_iterations = 200;

[M N] = size(Phi);
block_size = sqrt(N);

TOL = 0.0001;
D_prev = 0;

x = Phi' * y;

for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
      lambda, contourlet);

  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    break;
  end

  D_prev = D;
end

[x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
    lambda, contourlet, 'final');

reconstructed_image = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distict');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, block_size, num_rows, num_cols, ...
    lambda, contourlet, last)

x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct'); 
x_hat = wiener2(x, [3, 3]);
x_hat = im2col(x_hat, [block_size block_size], 'distinct');

x_hat = x_hat + Phi' * (y - Phi * x_hat);

x1 = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

x_check = pdfbdec(x1, contourlet.pfilter, contourlet.dfilter, ...
    contourlet.nlevels);
if (nargin == 9)
  start_level = length(x_check);
else
  start_level = length(x_check) - 2;
end
x_check = SPLBivariateShrinkage(x_check, start_level, lambda);
x_bar = double(pdfbrec(x_check, contourlet.pfilter, contourlet.dfilter));
x_bar = im2col(x_bar, [block_size block_size], 'distinct');

x = x_bar + Phi' * (y - Phi * x_bar);

x2 = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');
  
D = RMS(x1, x2);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x_check = SPLBivariateShrinkage(x_check, start_level, lambda)

windowsize  = 3;
windowfilt = ones(1, windowsize)/windowsize;

last_level = length(x_check);
num_last_level = length(x_check{last_level});
tmp = cell2mat(x_check{last_level}(1:floor(num_last_level/2)));
Nsig = median(abs(tmp(:)))/0.6745;

for i = last_level:-1:start_level
  num_cells_child = length(x_check{i});
  num_cells_parent = length(x_check{i-1});
  
  if num_cells_child ~= num_cells_parent
    for j = 1:floor(num_cells_parent/2)
      x_child = x_check{i}{j};
      x_parent = x_check{i-1}{j};
      
      x_child2 = x_check{i}{floor(num_cells_child/2)+j};
      x_parent2 = x_check{i-1}{floor(num_cells_parent/2)+j};
      
      x_parent1 = expand_col(x_parent);
      x_parent2 = expand_row(x_parent2);
      Wsig = conv2(windowfilt, windowfilt, (x_child).^2, 'same');
      Ssig = sqrt(max(Wsig-Nsig.^2, eps));
      
      T = sqrt(3)*Nsig^2./Ssig;
      x_child = bishrink(x_child, x_parent1, T*lambda);
      
      Wsig = conv2(windowfilt, windowfilt, (x_child2).^2, 'same');
      Ssig = sqrt(max(Wsig-Nsig.^2, eps));
      
      T = sqrt(3)*Nsig^2./Ssig;
      x_child2 = bishrink(x_child2, x_parent2, T*lambda);
      x_check{i}{j} = x_child;
      x_check{i}{floor(num_cells_child/2)+j} = x_child2;
    end
  else
    for j = 1:num_cells_parent
      x_child = x_check{i}{j};
      x_parent = x_check{i-1}{j};
      
      x_parent = expand(x_parent);
      Wsig = conv2(windowfilt, windowfilt, (x_child).^2, 'same');
      Ssig = sqrt(max(Wsig-Nsig.^2, eps));
      
      T = sqrt(3)*Nsig^2./Ssig;
      x_child = bishrink(x_child, x_parent, T*lambda);
      x_check{i}{j} = x_child;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = expand_col(x)

[N, M] = size(x);
M = M*2;

y = zeros(N, M);
y(1:N, 1:2:M) = x;
y(1:N, 2:2:M) = x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = expand_row(x)

[N, M] = size(x);
N = N*2;

y = zeros(N, M);
y(1:2:N, 1:M) = x;
y(2:2:N, 1:M) = x;
