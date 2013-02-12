% bpdq_proj_lpball_mex : Projection onto Lp ball via Newton's method
%
% [x_out,n] = bpdq_proj_lpball_mex(y_in,c,r,p)
%
% This function computes the projection of y_in onto the lp ball, ie
% solves 
%
% x_out = argmin ||x-y_in|| s.t. ||x-c||_p <= r
%
% Algorithm proceeds by multivariate Newton's method applied to
% find the zeros of the N+1 dimensional nonlinear equation F(z)=0
% arising from lagrange multiplier condition for the constrained 
% optimization problem.
%
% Algorithm should not be used for p<2 as Jacobian entries contain
% z_k^(p-2) ... 
%
% Inputs :
% y_in - input vector, to be projected
% c - center of Lp ball
% r - radius of Lp ball
% p - exponent for Lp ball
%
% Outputs:
% x_out - computed projection
% n - number of iterations
%
% Routine is implemented as a mex file. 
% For a more detailed explanation of the algorithm, see also the 
% matlab implementation bpdq_proj_lpball.m
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
