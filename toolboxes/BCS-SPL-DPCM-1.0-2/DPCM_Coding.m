% 
% function [y_q rate] = DPCM_Coding(y, quantizer_bitdepth, num_rows, num_cols)
% 
%   This function performs DPCM coding on the block measurements. For 
%   simplicity, this function assumes that coding and decoding occurs
%   at the same place. Therefore, the final output is recovered measurement.
%   Unlike BCS-SPL and MH-BCS-SPL, MS-BCS-SPL has different data structure 
%   and shows high correlation only on the baseband. Therefore, the procedure
%   is divided by two and handles the data accordingly.
%
% 
%   Inputs
%   y: block-by-block sampled measurements
%   quantizer_bitdepth: bit width
%   num_rows: number of rows of the original image
%   num_cols: number of columns of the original image
%
%   Outputs
%   y_q: quantized measurements
%   rate: rate information reprensted as entropy

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

function [y_q rate] = DPCM_Coding(y, quantizer_bitdepth, num_rows,num_cols)

if iscell(y) % quantization for cell form of data (MS-BCS-SPL)
  y_max = max(y{1}(:));   y_min = min(y{1}(:));

  % q: stepsize
  q = (y_max - y_min)/2^quantizer_bitdepth;
  
  yq = cell(1,length(y)); y_q = cell(1,length(y));
  
  yq{1}= DPCM_Encode(y{1},q(1));
  % y_rate is to measure entropy
  y_rate = yq{1}(:);
  y_q{1} = DPCM_Decode(yq{1},q(1));
  
  for i = 2:length(y)
    for subband = 1:3
      yq{i}{subband} = round(y{i}{subband}/q);

      y_rate = [y_rate; yq{i}{subband}(:)];      
      
      y_q{i}{subband} = yq{i}{subband}*q;
    end
  end
  total_pixels = num_rows*num_cols;
  rate = Measurement_Entropy(y_rate(:),total_pixels);

else % quantization for vector form or matrix form of data
  y_max = max(y(:));  y_min = min(y(:));
  q = (y_max - y_min)/2^quantizer_bitdepth;
  
  yq = DPCM_Encode(y,q);
  total_pixels = num_rows*num_cols;
  rate = Measurement_Entropy(yq(:),total_pixels);
  % yq can be coded using entropy coder, i.e. arithmetic coder. 
  
  y_q = DPCM_Decode(yq,q);
  
end


function i = DPCM_Encode(y,q)

i = zeros(size(y));
y_hat = zeros(size(y,1),1);

for j = 1:size(y,2)
    d = y(:,j) - y_hat;
    i(:,j) = round(d/q);
    d_hat = i(:,j)*q;
    y_hat = y_hat + d_hat;
end


function y_hat = DPCM_Decode(i,q)

y_hat = zeros(size(i));
y_hat_prev = zeros(size(i,1),1);
for j = 1:size(i,2)
    d_hat = i(:,j)*q;
    y_hat(:,j) = y_hat_prev + d_hat;
    y_hat_prev = y_hat(:,j); 
end
