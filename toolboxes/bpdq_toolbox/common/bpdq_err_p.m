% bpdq_err_p : Compute appropriate epsilon for BPDQ program
%
% epsilon=bpdq_err_p(p,alpha,M)
%
% ref : Jacques et al 2009, section III-C
%
% Inputs:
% p - BPDQ moment
% alpha - quantization bin width
% M  - # of measurements
%
% Outputs:
% epsilon - epsilon large enough so that original signal is feasible
% solution of BPDQ fidelity constraint, with high probability
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function epsilon=bpdq_err_p(p,alpha,M)
  kappa=2;
  epsilon = (alpha/2)*(M/(p+1))^(1/p)*(1+kappa*(p+1)/M^(1/2))^(1/p);
  
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
