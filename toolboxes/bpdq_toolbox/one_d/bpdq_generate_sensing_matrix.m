% bpdq_generate_sensing_matrix : Compute pseudorandom sensing matrix
%
% A=bpdq_generate_sensing_matrix(M,N,seed)
%
% Generates M x N Gaussian random matrix
%
% Inputs :
% M,N - dimensions of matrix
% seed - (optional) seed for random number generator
%
% Outputs :
% A - M x N matrix with zero mean, unit variance i.i.d. Gaussian entries.
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function A=bpdq_generate_sensing_matrix(M,N,seed)
  if (nargin==3)
    randn('seed',seed);
  end
  A = randn(M,N);
  
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
