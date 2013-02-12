%
% function run_experiment_msbcsspl_dpcm()
%
%   This function runs the MS-BCS-SPL-DPCM algorithm.

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

function run_experiment_msbcsspl_dpcm()

addpath(genpath('../BCS-SPL-1.5-1'));
addpath('../WaveletSoftware');
addpath(genpath('../MS-BCS-SPL-1.2-1'));
addpath('../waveletcdf97');

filename = 'lenna';
% filename = 'barbara';
% filename = 'goldhill';
% filename = 'peppers';
% filename = 'clown';
% filename = 'man';

original_image = double(imread([filename '.pgm']));
[num_rows num_cols] = size(original_image);
x = original_image;

block_sizes = [2 2 2];
num_levels = length(block_sizes);

% load bits and subrates used in the paper 
algorithm = 'msbcsspl';
type = 'sq';
filename_parameters_SQ = ['./Parameters/' filename '_' algorithm '_' ...
      type '_parameters.mat' ];
load(filename_parameters_SQ);

num_trials = index;
rate_SQ = zeros(1, num_trials);
psnr_SQ = zeros(1, num_trials);
for trial = 1:num_trials
  quantizer_bitdepth = bits(trial);
  subrate = subrates(trial);

  subrates_ms = MS_BCS_SPL_GetSubrates(subrate, num_levels);
  Phi = MS_BCS_SPL_GenerateProjection(block_sizes, subrates_ms, ...
      ['projections.' num2str(num_levels) '.mat'], 'dct');
  
  % random projection
  y = MS_BCS_SPL_Encoder(x, Phi);
  
  % uniform scalar quantization
  [y_SQ rate_SQ(trial)] = SQ_Coding(y, quantizer_bitdepth, num_rows, num_cols);
  
  % reconstruction
  disp('Decoding SQ...');
  tic
  x_hat_SQ = MS_BCS_SPL_DDWT_Decoder(y_SQ, Phi, subrates_ms, num_rows, ...
      num_cols, num_levels);
  toc
  psnr_SQ(trial) = PSNR(original_image, x_hat_SQ);
  disp(['Rate = ' num2str(rate_SQ(trial), '%0.4f') ...
	' (bpp), PSNR = ' num2str(psnr_SQ(trial), '%0.2f') ' (dB)']);
end

% load bits and subrates used in the paper 
algorithm = 'msbcsspl';
type = 'dpcm';
filename_parameters_SQ = ['./Parameters/' filename '_' algorithm '_' ...
      type '_parameters.mat' ];
load(filename_parameters_SQ);

num_trials = index;
rate_DPCM = zeros(1, num_trials);
psnr_DPCM = zeros(1, num_trials);
for trial = 1:num_trials
  quantizer_bitdepth = bits(trial);
  subrate = subrates(trial);

  subrates_ms = MS_BCS_SPL_GetSubrates(subrate, num_levels);
  Phi = MS_BCS_SPL_GenerateProjection(block_sizes, subrates_ms, ...
      ['projections.' num2str(num_levels) '.mat'], 'dct');
  
  % random projection
  y = MS_BCS_SPL_Encoder(x, Phi);
  
  % DPCM + uniform scalar quantization
  [y_DPCM rate_DPCM(trial)] = ...
      DPCM_Coding(y, quantizer_bitdepth, num_rows, num_cols);
  
  % reconstruction
  disp('Decoding DPCM...');
  tic
  x_hat_DPCM = MS_BCS_SPL_DDWT_Decoder(y_DPCM, Phi, subrates_ms, ...
      num_rows, num_cols, num_levels);
  toc
  
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

legend('MS-BCS-SPL+DPCM', 'MS-BCS-SPL+SQ', 'Location', 'SouthEast');

print('-deps', [filename '_' algorithm '.eps']);
