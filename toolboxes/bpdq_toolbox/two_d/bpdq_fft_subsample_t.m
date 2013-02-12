% bpdq_fft_subsample_t : Transpose of operation bpdq_fft_subsample
%
% r=bpdq_fft_subsample_t(x,ind,dim)
% Computes transpose of fft_subsample
%
% Inputs:
% x - input signal (complex vector same size as ind)
% ind - indices of desired Fourier coefficients
% dim - image size
% 
% Outputs:
% r - returned transpose, will be vector of size prod(dim) x 1
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function r=bpdq_fft_subsample_t(x,ind,dim)
  N=length(x);
  rdft=zeros(dim);
  rdft(ind)=x;
  r=ifft2(ifftshift(rdft));
  r=r*prod(dim); % fft normalization 
  r=r(:);
  
% The BPDQ Toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%  
% The BPDQ Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License
% along with The BPDQ Toolbox.  If not, see <http://www.gnu.org/licenses/>.
