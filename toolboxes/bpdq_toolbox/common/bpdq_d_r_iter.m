% bpdq_d_r_iter : Douglas-Rachford Iteration
%
% [xstar,relerr_save]=bpdq_d_r_iter(proxf1,proxf2,gamma,lambda,x0,varargin)
%
% Implements Douglas-Rachford Splitting to solve
% xstar = argmin_x [ f1(x) + f2(x) ]
%
% Inputs:
% proxf1, proxf2 - function handles computing prox(gamma*f1), prox(gamma*f2)
% i.e., proxf1(y,gamma) solves argmin_x [ (1/2) |y-x|^2 + gamma * f(x) ]
% gamma, lambda -  Douglas-Rachford algorithm parameters
% x0 - starting point
%
% Selectable control parameters:
% dr_maxiter - maximum # of iterations
% dr_relerr_thresh - relative error threshold to terminate iteration
% dr_verbosity - 0 to supress output, 1 to output relative error at each 
% iteration
%
% Outputs :
% xstar - returned solution
% relerr_save - record of relative error of iterates; 
% relerr_save(n) = ||x_n-x_(n-1)||/||x_n||
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function [xstar,relerr_save]=bpdq_d_r_iter(proxf1,proxf2,gamma,lambda,...
                                           x0,varargin)
  control_params={'dr_maxiter',500,...
                  'dr_relerr_thresh',1e-6,...
                  'dr_verbosity',0};
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);

  xp_old=zeros(size(x0));
  xp_new=zeros(size(x0));
  x_new=zeros(size(x0));
  x_old=x0;
  for n=1:dr_maxiter
    xp_old=xp_new; 
    if (n>1)
      x_old=x_new;
    end
    
    xp_new=proxf1(x_old,gamma);
    x_new = x_old + lambda*(proxf2(2*xp_new-x_old,gamma)-xp_new);
    relerr = norm(xp_old-xp_new)/norm(xp_new);
    relerr_save(n)=relerr;
    if (dr_verbosity>=1)
      fprintf('DR iter %g, relerr = %g\n',n,relerr);
    end
    if (relerr<dr_relerr_thresh) && (n>3)
      break;
    end      
    %keyboard
  end
  xstar=proxf1(x_new,gamma);
  
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
