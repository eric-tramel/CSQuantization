% bpdq_proj_lpball : Projection onto Lp ball via Newton's method
%
% function [x_out,n] = bpdq_proj_lpball(y_in,c,r,p)
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
% This is a matlab implementation, if possible the compiled c++
% implementation bpdq_proj_lpball_mex should be used instead.
% This matlab file may be considered a more easily readable reference
% implementation of the algorithm, and is included largely for 
% pedagogic purposes
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function [x_out,n] = bpdq_proj_lpball(y_in,c,r,p)
  tol=1e-15;
  maxiter=50;
  verbose=0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check first if input is inside specified ball
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if norm(y_in - c,p)<r
    x_out=y_in;
    n=0;
    return
  end
  
  if (p == Inf)
    x_out = c + sign(y_in - c).*min(abs(y_in - c), r);
    n=0;
    return
  end
  
  if p<2
    error('Don''t try to use this code for p<2');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % preprocessing to transform to c=0; r=1; case
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y=(y_in-c)/r;
  % preprocessing to remove signs 
  % (work in orthant where everything is non-negative)
  ysign=sign(y);
  y=y.*ysign;
  
  % now solve for xstar for c=0, r=1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initialization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  N=numel(y);
  F=zeros(N+1,1);
  z=zeros(N+1,1);
  % xn1 : radial projection onto Lp ball
  xn1=y/norm(y,p); 
  % xn2 : L\infty projection followed by radial projection onto Lp ball
  xn2=sign(y).*min(abs(y),1);
  xn2=xn2/norm(xn2,p);
  if (norm(xn1-y) < norm(xn2-y))
    x0=xn1;
  else
    x0=xn2;
  end  
  z(1:N)=x0;
  % pick least squares for lambda solving 
  % lambda * x_i^(p-1) = y_i-x_i
  % as best it can 
  lambda = x0.^(p-1) \ (y-x0);
  z(N+1)=lambda;
  
  % Indices for constructing sparse Jacobian
  diag_i=1:N;
  diag_j=1:N;
  rightcol_i=1:N;
  rightcol_j=(N+1)*ones(1,N);
  bottomrow_i=(N+1)*ones(1,N);
  bottomrow_j=1:N;
  
  all_i=[diag_i(:);rightcol_i(:);bottomrow_i(:)];
  all_j=[diag_j(:);rightcol_j(:);bottomrow_j(:)];
  
  normz0=norm(z); % norm of input, useful for computing relative tolerance
  normF=tol+1;
  n=1;
  xn=z(1:N);
  % begin iteration

  while normF/normz0 > tol;
    % Build residual F
    zpm1=z(1:N).^(p-1); % z to the p minus 1
    F(1:N)=z(1:N)+z(N+1)*zpm1-y;
    % divide by p to yield symmetric Jacobian ( ref my notes Nov 19)

    F(N+1)= (sum(z(1:N).^p)-1)/p;
    normF=norm(F);

    % Build Jacobian matrix J
    d=1+z(N+1)*(p-1)*z(1:N).^(p-2); % diagonal part
        
    %%
    % This code computes the equivalent of inv(J)*-F using
    % the block inverse formula for J = [D,b ; b',0]
    % where D is diagonal.
    %%
    btwid=(zpm1./d); % inv(A)*b
    mu=zpm1'*btwid; % b^T*inv(A)*b
    chi=btwid'*(-F(1:N));
    dz=[ [-F(1:N)./d+(-F(N+1)-chi)/mu*btwid]; (chi-(-F(N+1)) )/mu ];      
    
    z=z+dz;
    
    if verbose
      oxn=xn;
      xn=z(1:N);
      dxn=norm(xn - oxn)/norm(xn);
      cosn = dot(lpnormalvector(xn,p), (y- xn)/norm(y - xn));
      fprintf('%i: |y-xn|_2=%e, |xn+1-xn|_2/|xn+1|=%e, cosn=%e, nF=%e\n',...
              n,norm(y - xn), dxn,cosn,normF);
    end
    n=n+1;
    if n>maxiter
      error('maximum # of iterations exceeded in proj_lpball_newton');
    end
  end
  
  xstar=z(1:N);
  % postprocessing to restore signs
  xstar=xstar.*ysign;
  
  % postprocessing to transform back to c!=0; r!=1 case
  x_out=r*xstar+c;
 
  
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
