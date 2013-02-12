% prox_tv_mex : compute prox of TV operator
%
% xstar = prox_tv_mex(y,lambda,...)
%
% Solves xstar= min 0.5||y-x||_2 + lambda ||x||_TV
%
% The problem is solved thanks to the forward-backward method described 
% in: 
% Multiplicative Noise Removal Using L1 Fidelity on Frame Coefficients,
%  Durand S., Fadili J. and Nikolova M.
%  arXiv:0812.1697v1
%
% Inputs : 
% y - input image
% lambda - regularization parameter
%
% Optional parameters (must be passed in 'name',value pairs)
% beta (0.249) - parameter for forward-backward splitting algorithm
% min_rel_obj (0.001) - 
% algorithm will terminate if minimum relative objective is below this value
% it_max (500) - maximum # of iterations
% it_min (1)  - minimum # of iterations
% verbose (0) - whether to print output at each iteration
%
% Outputs:
% xstar - returned solution
%
% Routine is implemented as mex file
% 
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

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
