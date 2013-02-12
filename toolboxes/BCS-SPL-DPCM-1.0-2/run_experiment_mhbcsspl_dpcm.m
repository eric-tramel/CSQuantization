% 
% function run_experiment_mhbcsspl_dpcm()
% 
%   This function runs the MH-BCS-SPL-DPCM algorithm.
%

%
% DPCM for BCS-SPL: Differential Pulse-Code Modulation for 
% Block Compressed Sensing - Smooth Projected Landweber
% Copyright (C) 2012  James E. Fowler
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

%
%   Originally written by SungKwang Mun, Mississippi State University
%


function run_experiment_mhbcsspl_dpcm()

addpath(genpath('../BCS-SPL-1.5-1'));
addpath(genpath('../MH-BCS-SPL-1.0-1'));
addpath(genpath('../WaveletSoftware'));
addpath(genpath('../SSIM'));

filename = 'lenna';
% filename = 'barbara';
% filename = 'goldhill';
% filename = 'peppers';
% filename = 'clown';
% filename = 'man';

original_image = double(imread([filename '.pgm']));
[num_rows num_cols] = size(original_image);
x = original_image;

block_size = 16;
projection_matrix_file = ['projections.' num2str(block_size) '.mat'];

% load bits and subrates used in the paper 
algorithm = 'mhbcsspl';
type = 'sq';
filename_parameters_SQ = ['./Parameters/' filename '_' algorithm '_' ...
      type '_parameters.mat' ];
load(filename_parameters_SQ);

num_trials = length(subrates);
rate_SQ = zeros(1, num_trials);
psnr_SQ = zeros(1, num_trials);
for trial = 1:num_trials
  quantizer_bitdepth = bits(trial);
  subrate = subrates(trial);
  
  % random projection
  Phi = BCS_SPL_GenerateProjection(block_size, subrate, ...
      projection_matrix_file);
  y = BCS_SPL_Encoder(x, Phi);
  
  % uniform scalar quantization
  [y_SQ rate_SQ(trial)] = SQ_Coding(y, quantizer_bitdepth, num_rows, num_cols);
  
  disp('Decoding SQ...');
  x_hat_SQ = MH_BCS_SPL_Decoder(y_SQ, Phi, subrate, num_rows, num_cols);
  psnr_SQ(trial) = PSNR(original_image, x_hat_SQ);
  disp(['Rate = ' num2str(rate_SQ(trial), '%0.4f') ...
	' (bpp), PSNR = ' num2str(psnr_SQ(trial), '%0.2f') ' (dB)']);
end

% load bits and subrates used in the paper 
algorithm = 'mhbcsspl';
type = 'dpcm';
filename_parameters_SQ = ['./Parameters/' filename '_' algorithm '_' ...
      type '_parameters.mat' ];
load(filename_parameters_SQ);

num_trials = length(subrates);
rate_DPCM = zeros(1, num_trials);
psnr_DPCM = zeros(1, num_trials);
for trial = 1:num_trials
  quantizer_bitdepth = bits(trial);
  subrate = subrates(trial);
  
  % random projection
  Phi = BCS_SPL_GenerateProjection(block_size, subrate, ...
      projection_matrix_file);
  y = BCS_SPL_Encoder(x, Phi);
  
  % DPCM + uniform scalar quantization
  [y_DPCM rate_DPCM(trial)] = ...
      DPCM_Coding(y, quantizer_bitdepth, num_rows, num_cols);
  
  % reconstruction
  disp('Decoding DPCM...');
  x_hat_DPCM = MH_BCS_SPL_Decoder(y_DPCM, Phi, subrate, num_rows, num_cols);
  psnr_DPCM(trial) = PSNR(original_image, x_hat_DPCM);
  disp(['Rate = ' num2str(rate_DPCM(trial), '%0.4f') ...
	' (bpp), PSNR = ' num2str(psnr_DPCM(trial), '%0.2f') ' (dB)']);
end

save([filename '_' algorithm '_results.mat'], ...
    'rate_SQ', 'psnr_SQ', ...
    'rate_DPCM', 'psnr_DPCM');

figure(1);
clf;
plot(rate_DPCM, psnr_DPCM, 'k-x', 'LineWidth', 2);
hold on
plot(rate_SQ, psnr_SQ, 'k--x', 'LineWidth', 2);
grid on
title(filename);
xlabel('Bitrate (bpp)');
ylabel('PSNR (dB)');

legend('MH-BCS-SPL+DPCM', 'MH-BCS-SPL+SQ', 'Location', 'SouthEast');

print('-deps', [filename '_' algorithm '.eps']);
