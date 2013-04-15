function [S B] = subrate_bit_LUT(bitrate)
%SUBRATE_BIT_LUT returns subrate/bit-depth pair for a target bitrate.
%			
% function [S B] = subrate_bit_LUT(bitrate)
%
%  This function will return a subrate (S) and bit-depth (B) for a
%  given floating-point bitrate target. BE ADVISED: This LUT has defined
%  values within the range [0.1,1.5]bpp. If you are using values outside
%  that range, either update the LUT in this funciton or use AT YOUR OWN
%  PERIL.
%
%
% Written by: Eric W. Tramel, Ph.D.
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

subrates = [0.05 0.08 0.11 0.11 0.14 0.17 0.2 0.25 0.33 0.37 0.4 0.44 0.47 0.51 0.54];
bits = [4 5 5 6 6 6 6 6 6 6 6 6 6 6 6]; 
bitrates = linspace(0.1,1.5,15);

S = interp1(bitrates,subrates,bitrate,'linear','extrap');
B = max(round(interp1(bitrates,bits,bitrate,'linear','extrap')),1);