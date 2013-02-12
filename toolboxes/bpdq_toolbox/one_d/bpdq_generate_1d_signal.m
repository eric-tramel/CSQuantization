% bpdq_generate_1d_signal : Generate random sparse signal
%
% x = bpdq_generate_1d_signal(N,K,seed)
% Generates N length K sparse random signal 
% (K nonzero entries are iid zero mean Gaussian  with unit variance)
%
% Inputs :
% N - Signal length
% K - sparsity
% seed - (optional) seed for random number generator
%
% Outputs:
% x- returned random signal
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function x = bpdq_generate_1d_signal(N,K,seed)
  if (nargin==3)
    rand('seed',seed);
    randn('seed',seed);  
  end
  x = zeros(N,1);
  pos = randperm(N);
  x(pos(1:K)) = randn(K,1);
  
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
