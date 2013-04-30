% 
% function H = Measurement_Entropy(y, total_pixels)
% 
%   This function performs calculation of the entropy of the input 
%   frequency of the measurements y. Because of the dimensional reduction 
%   property of CS, we take into account the number of total pixels in
%   the image.
%   
%   Inputs
%   y: (quantized) measurements
%   total_pixels: total number of pixels of the image to be reconstructed
%
%   Outputs
%   H: entropy
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

function H = Measurement_Entropy(y, total_pixels)

% assumed that y is vector.
y = y(:);

nbin = unique(y);
frequency = hist(y,nbin);

p = frequency/sum(frequency);
p(p==0) = 1;
H = sum(-p.*log2(p));

% total pixel is used to calculate the real entropy
H = H*numel(y)/total_pixels;

  
