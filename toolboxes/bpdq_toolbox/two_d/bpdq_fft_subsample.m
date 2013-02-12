% bpdq_fft_subsample : Linear operation of Fourier transform and subsampling
%
% r=bpdq_fft_subsample(x,ind,dim)
%
% Implements the linear operation corresponding to subsampling 2d Fourier
% transform at locations specified by ind.
%
% Output is vectorized, as complex values
%
% Inputs:
% x - input signal (in spatial domain)
% ind - indices of desired Fourier coefficients (with FT shifted by fftshift
% so that DC coefficient is in center of image)
% dim - image size
%
% Outputs:
% r - returned sampled Fourier values, returned as complex vector
% the same size as input ind
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function r=bpdq_fft_subsample(x,ind,dim)
  xdft=fftshift(fft2(reshape(x,dim)));
  r=xdft(ind);

  
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
