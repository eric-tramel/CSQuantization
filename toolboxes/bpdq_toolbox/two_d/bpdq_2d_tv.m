% bpdq_2d_tv : Implementation of BPDQ-TV for 2d reconstruction 
%
% [xstar,D] = bpdq_2d_tv(yq,A,At,dim,epsilon,p,varargin)
%
% Implements reconstruction from quantized Fourier measurements
%
% xstar  = argmin_x ||x||_TV s.t. ||yq-A*x||_p <= epsilon
%
% Inputs:
% yq - quantized measurements
% A - function handle implementing sensing matrix
% At - function handle implementing transpose of sensing matrix
% dim - 2d dimensions of image 
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function [xstar,D] = bpdq_2d_tv(yq,A,At,dim,epsilon,p,varargin)
  control_params={'dr_verbosity',0,...
                  'dr_maxiter',50,...
                  'dr_gamma',.1,...
                  'dr_lambda',1};
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);
  
  % unpack complex vector into vector of reals
  c2r = @(x) [real(x(:));imag(x(:))];
  r2c = @(x) x(1:end/2)+sqrt(-1)*x(end/2+1:end);
  
  % Assuming A*A' = nu * I ; compute nu
  x=rand(size(yq));
  nu=(At(x))'*At(x)/(x'*x);
  
  % define prox operators
  % f1 is modified lpball. As A is tight frame, fca type iteration
  % is not needed, can use direct formula to get prox(f(A(x))
  prox_lp =@(x) r2c(bpdq_proj_lpball_mex( c2r(x),c2r(yq),epsilon,p));
  proxf1 =@(x,tgamma) x+ (1/nu)*At( prox_lp(A(x)) - A(x));  
  
  % f2 is TV norm. Note : the need to reshape back to 2-d array to compute
  % prox of TV norm is why dim must given
  tvprox_opts={'min_rel_obj',5e-4,'it_max',500,'verbose',0};
  tv_prox=@(x,tgamma)vec(bpdq_prox_tv_mex(real(reshape(x,dim)),tgamma,...
                                          tvprox_opts{:}));
  proxf2=tv_prox;
  
  % init point
  % x0=real( vec( samp_t(yq))/prod(dim) );
  x0=At(yq)/prod(dim);
  [xstar,dr_relerr_save]=bpdq_d_r_iter(proxf1,proxf2,dr_gamma,dr_lambda,x0,...
                                       'dr_verbosity',dr_verbosity,...
                                       'dr_maxiter',dr_maxiter);
  D.dr_relerr_save=dr_relerr_save;
  
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
