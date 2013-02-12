% bpdq_prox_f_comp_A : compute prox of g(x) = f(A(x))
%
% function [r,relerr_save] = bpdq_prox_f_comp_A(x,proxcf,opA,opAt,c,varargin)
%
% Prox is computed by iteration (lemma 4 of Jacques et al)
%
% Inputs:
% x - input independent variable
% proxcf - function handle for computing prox_(c*f)
% opA, opAt - function handles computing multiplication by A and A transpose
% c - constant, should be greater than half of operator norm of A
%
% Selectable control parameters :
% fca_maxiter - maximum # of iterations
% fca_relerr_thresh - relative error threshold to terminate iteration
% fca_verbosity - 0 to supress output, 1 to output # of iterations until
% threshold reached, 2 to output relative error at each iteration.
%
% Outputs : 
% r - prox of g(x)=f(A(x))
% relerr_save - record of relative error of iterates.
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function [r,relerr_save] = bpdq_prox_f_comp_A(x,proxcf,opA,opAt,c,varargin)
  control_params={'fca_maxiter',500,...
                  'fca_relerr_thresh',1e-6,...
                  'fca_verbosity',0};
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);  

  u_old = opA(x);
  p_old=x;
  for n=1:fca_maxiter 
    p_new = x-opAt(u_old);
    tmp1=c*u_old + opA(p_new);
    u_new = (1/c)*(tmp1-proxcf(tmp1,c));

    relerr=norm(p_new-p_old)/norm(p_new);
    relerr_save(n)=relerr;

    u_old=u_new;
    p_old=p_new;
    if (fca_verbosity>=2)
      fprintf('fca iteration n %g, relerr %g\n',n,relerr);
    end
    if (relerr<fca_relerr_thresh)
      if (fca_verbosity>=1)
        fprintf('fca relerr %g reached at %g iterations\n',...
                relerr,n);
      end
      break;
    end    
  end
  r=x-opAt(u_new);
  
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
