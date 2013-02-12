% bpdq_quantize : Uniform scalar quantization
%
% xq = bpdq_quantize(x,alpha)
%
% Performs uniform quantization with bin width alpha
% aligned so that bin edges are 
% ... -2*alpha, -alpha , 0 , alpha ,2*alpha ,...
% 
% Scalar quantization is applied independently to each component of 
% input vector.
%
% Quantized values take center of bin to which they belong
% e.g. so that Qd(0.1,1)=0.5
%
% Input may be complex, in which case quantization is applied to
% real and imaginary parts independently
%
% Inputs :
% x - vector input
% alpha - scalar quantization bin width
%
% Outputs :
% xq - quantized output
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function xq = bpdq_quantize(x,alpha)
  if isreal(x)
    xq=alpha*floor(x/alpha) + alpha/2;
  else
    xq=bpdq_quantize(real(x),alpha)+i*bpdq_quantize(imag(x),alpha);
  end

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
