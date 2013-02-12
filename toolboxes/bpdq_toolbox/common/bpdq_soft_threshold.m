% bpdq_soft_threshold : Soft thresholding operator
%
% x_t = bpdq_soft_threshold(x,tgamma)
%
% Applies soft thresholding  to each component of x
%
% Inputs:
% x - input signal
% tgamma - threshold
%
% Outputs:
% x_t - soft thresholded result
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function x_t = bpdq_soft_threshold(x,tgamma)
  tmp=abs(x)-tgamma;
  x_t = sign(x).*tmp.*(tmp>0);
  
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
