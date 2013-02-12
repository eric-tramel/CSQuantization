% 
% function psnr = run_experiment_ddwt()
% 
%   This function runs the experiments for BCS-SPL-DDWT in Table 1 of
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


function psnr = run_experiment_ddwt()

filename = 'lenna';
% filename = 'barbara';
% filename = 'goldhill';

num_trials = 5;
subrate = 0.3;
block_size = 32;
num_levels = 6;

addpath('Images');
addpath('../WaveletSoftware');

original_filename = [filename '.pgm'];
original_image = double(imread(original_filename));

[num_rows num_cols] = size(original_image);

for trial = 1:num_trials
  projection_matrix_file = ['projections.' num2str(block_size) '.' ...
        num2str(trial) '.mat'];
  
  Phi = BCS_SPL_GenerateProjection(block_size, subrate, projection_matrix_file);
  
  y = BCS_SPL_Encoder(original_image, Phi);
  
  reconstructed_image = BCS_SPL_DDWT_Decoder(y, Phi, num_rows, ...
      num_cols, num_levels);
  
  psnr(trial) = PSNR(original_image, reconstructed_image);
end

psnr = mean(psnr);

disp(['BCS-SPL-DDWT PSNR = ' num2str(psnr) ' dB']);

