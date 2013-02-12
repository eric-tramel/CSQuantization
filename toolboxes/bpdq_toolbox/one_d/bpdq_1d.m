% bpdq_1d : Main implementation of BPDQ algorithm
%
% x_bpdq = bpdq_1d(yq,A,epsilon,p,varargin);
%
% This solves the BPDQ program
% x_bpdq = argmin_x  ||x||_1 s.t. ||yq - A*x||_p <= epsilon
% using the  Douglas-Rachford algorithm.
%
% This is the implementation of the method described in
% "Dequantizing Compressed Sensing : When Oversampling and Non-Gaussian
% Constraints Combine"
% by Laurent Jacques, David Hammond, M. Jalal Fadili
% ( submitted to IEEE Transactions on Signal Processing)
%
% Inputs:
% yq - quantized measurements
% A - sensing matrix
% epsilon - lpball radius, chosen so that original signal is 
%   feasible with high probability
% p - BPDQ moment
%
% Selectable control parameters :
% dr_gamma, dr_lambda - DR algorithm parameters
% dr_maxiter - maximum number of DR iterations
% dr_relerr_thresh - relatize error threshold to terminate DR iteration
% dr_verbosity - 1 to enable output, 0 to supress
% fca_maxiter - max number of iterations for fca iteration
% (f composed with A - iteration to compute prox(f(A(x))) : see 
% bpdq_f_comp_A.m for details
% fca_relerr_thresh - relative error threshold to terminate fca iteration
% fca_verbosity - 1 to enable output, 0 to supress
%
% Outputs:
% x_bpdq - solution to BPDQ program
% D - debugging structure containing field dr_relerr_save, for
% verifying empirical convergence of algorithm
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function [x_bpdq,D] = bpdq_1d(yq,A,epsilon,p,varargin)
  control_params={'dr_gamma',.05,...
                  'dr_lambda',.5,...               
                  'dr_maxiter',500,...
                  'dr_relerr_thresh',1e-6,...
                  'dr_verbosity',0,...
                  'fca_maxiter',500,...
                  'fca_relerr_thresh',1e-6,...
                  'fca_verbosity',0};
  
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);
  
  % c is upper frame bound for frame corresponding to matrix A
  % (used for iteration for composition, ref
  % Jacques et. al. Lemma 4)
  c= max(svd(A))^2;
  
  % Sensing data
  [M,N] = size(A);
  opA = @(in) A*in;
  opAt = @(in) A'*in;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % define prox operators
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % f1 is ||x||_1, prox is soft thresholding
  prox1= @(u,tgamma) bpdq_soft_threshold(u,tgamma);
    
  % f2 is ||yq-A*x||_p, prox computed iteratively
  % with function composition and lp_ball projection
  % (tgamma may be discarded as f1(x) has only zero and infinite values)
  lpball_prox = @(u,tgamma) bpdq_proj_lpball_mex(u,yq,epsilon,p);
  prox2 = @(u,tgamma) bpdq_prox_f_comp_A(u,lpball_prox,opA,opAt,c,...
      'fca_maxiter',fca_maxiter,...
      'fca_relerr_thresh',fca_relerr_thresh,...
      'fca_verbosity',fca_verbosity);
  
  x0=zeros(N,1); % initialize at zero
  % call douglas rachford iteration with defined prox operators
  [x_bpdq,relerr_save]=bpdq_d_r_iter(prox1,prox2,dr_gamma,dr_lambda,x0,...
      'dr_maxiter',dr_maxiter,...
      'dr_relerr_thresh',dr_relerr_thresh,...
      'dr_verbosity',dr_verbosity);
  D.dr_relerr_save=relerr_save;

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
