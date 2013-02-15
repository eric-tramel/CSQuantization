function [theta,theta_debias,objective,times,debias_start,mses]= ...
    GPSR_BB(y,A,tau,varargin)
  
% 
%
% 
% This function provides a solver the convex problem 
% arg min_theta = 0.5*|| y - A theta ||_2^2 + tau ||theta||_1
% using the algorithm GPSR-BB, described in the paper
% "Gradient Projection for Sparse Reconstruction: Application
% to Compressed Sensing and Other Inverse Problems"
% by Mario A. T. Figueiredo, Robert D. Nowak, Stephen J. Wright
% submitted, 2007.
%
% Contacts: mario.figueiredo@lx.it.pt
%           nowak@ece.wisc.edu
%           swright@cs.wisc.edu
%
% Copyright (2007): Mario Figueiredo, Robert Nowak, Stephen Wright
%
% This code is in development stage; any comments or bug reports 
% are welcome.
%
%  ===== Required inputs =============
%
%  y: 1D vector or 2D array (image) of observations
%     
%  A: if y and theta are both 1D vectors, A can be a 
%     k*n (where k is the size of y and n the size of theta)
%     matrix or a handle to a function that computes
%     products of the form A*v, for some vector v.
%     In any other case (if y and/or theta are 2D arrays), 
%     A has to be passed as a handle to a function which computes 
%     products of the form A*x; another handle to a function 
%     AT which computes products of the form A'*x is also required 
%     in this case. The size of theta is determined as the size
%     of the result of applying AT.
%
%  tau: usually, a non-negative real parameter of the objective 
%       function (see above). It can also be an array, the same 
%       size as theta, with non-negative entries; in this  case,
%       the objective function weights differently each element 
%       of theta, that is, it becomes
%       0.5*|| y - A theta ||_2^2 + tau^T * abs(theta)
%
%  ===== Optional inputs =============
%
%  
%  'AT'    = function handle for the function that implements
%            the multiplication by the conjugate of A, when A
%            is a function handle. If A is an array, AT is ignored.
%
%  'StopCriterion' = type of stopping criterion to use
%                    0 = algorithm stops when the relative 
%                        change in the number of non-zero 
%                        components of the estimate falls 
%                        below 'ToleranceA'
%        any non-zero = algorithm stops when the relative 
%                       change in the objective function 
%                       falls below 'ToleranceA', for the
%                       first phase and 'ToleranceD' for 
%                       debiasing phase.
%                       Default = 0
%
%  'ToleranceA' = stopping threshold; Default = 0.01
% 
%  'Debias'     = debiasing option: 1 = yes, 0 = no.
%                 Default = 0.
%
%  'ToleranceD' = stopping threshold for the debiasing phase:
%                 Default = 0.0001.
%                 If no debiasing takes place, this parameter,
%                 if present, is ignored.
%
%  'MaxiterA' = maximum number of iterations allowed in the
%               main phase of the algorithm.
%               Default = 200
%
%  'MiniterA' = minimum number of iterations performed in the
%               main phase of the algorithm.
%               Default = 5
%
%  'MaxiterD' = maximum number of iterations allowed in the
%               debising phase of the algorithm.
%               Default = 200
%
%  'MiniterD' = minimum number of iterations to perform in the
%               debiasing phase of the algorithm.
%               Default = 5
%
%  'Initialization' must be one of {0,1,2,array}
%               0 -> Initialization at zero. 
%               1 -> Random initialization.
%               2 -> initialization with A'*y.
%           array -> initialization provided by the user.
%               Default = 0;
%             
%  'TrueTheta' = if the true underlying theta is passed in 
%                this argument, MSE plots are generated.
%
%  'AlphaMin' = the alphamin parameter of the BB method.
%               (See the paper for details)
%               Default = 1e-30;
%
%  'AlphaMax' = the alphamax parameter of the BB method.
%               (See the paper for details)
%               Default = 1e30;
%
% ===================================================  
% ============ Outputs ==============================
%   theta = solution of the main algorithm
%
%   theta_debias = solution after the debiasing phase;
%                  if no debiasing phase took place, this
%                  variable is empty, theta_debias = [].
%
%   objective = sequence of values of the objective function
%
%   times = CPU time after each iteration
%
%   debias_start = iteration number at which the debiasing 
%                  phase started. If no debiasing took place,
%                  this variable is returned as zero.
%
%   mses = sequence of MSE values, with respect to TrueTheta,
%          if it was given; if it was not given, mses is empty,
%          mses = [].
% ========================================================

% test for number of required parametres
if (nargin-length(varargin)) ~= 3
     error('Wrong number of required parameters');
end

% Set the defaults for the optional parameters
stop_objective = 0;
tolA = 0.01;
tolD = 0.0001;
debias = 0;
maxiter = 200;
maxiter_debias = 200;
miniter = 5;
miniter_debias = 5;
init = 0;
alphamin = 1e-30;
alphamax = 1e30;
compute_mse = 0;
AT = 0;

% Set the defaults for outputs that may not be computed
debias_start = 0;
theta_debias = [];
mses = [];

% Read the optional parameters
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch varargin{i}
     case 'StopCriterion'
       stop_objective = varargin{i+1};
     case 'ToleranceA'       
       tolA = varargin{i+1};
     case 'ToleranceD'
       tolD = varargin{i+1};
     case 'Debias'
       debias = varargin{i+1};
     case 'MaxiterA'
       maxiter = varargin{i+1};
     case 'MaxiterD'
       maxiter_debias = varargin{i+1};
     case 'MiniterA'
       miniter = varargin{i+1};
     case 'MiniterD'
       miniter_debias = varargin{i+1};
     case 'Initialization'
       if prod(size(varargin{i+1})) > 1   % we have an initial theta
          init = 33333;    % some flag to be used below
          theta = varargin{i+1};
       else 
          init = varargin{i+1};
       end
     case 'TrueTheta'
       compute_mse = 1;
       true = varargin{i+1};
     case 'AlphaMin'
       alphamin = varargin{i+1};
     case 'AlphaMax'
       alphamax = varargin{i+1};
     case 'AT'
       AT = varargin{i+1};
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%

% if A is a function handle, we have to check presence of AT,
if isa(A, 'function_handle') & ~isa(AT,'function_handle')
   error(['The function handle for transpose of A is missing']);
end 

% if A is a matrix, we find out dimensions of y and theta,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~isa(A, 'function_handle')
   AT = @(x) A'*x;
   A = @(x) A*x;
end
% from this point down, A and AT are always function handles.

% start the clock
t0 = cputime;

% Precompute A'*y since it'll be used a lot
Aty = AT(y);

% Initialization
switch init
    case 0   % initialize at zero, using AT to find the size of theta
       theta = AT(zeros(size(y)));
    case 1   % initialize randomly, using AT to find the size of theta
       theta = randn(size(AT(zeros(size(y)))));
    case 2   % initialize theta0 = A'*y
       theta = Aty; 
    case 33333
       % initial theta was given as a function argument; just check size
       if size(A(theta)) ~= size(y)
          error(['Size of initial theta is not compatible with A']); 
       end
    otherwise
       error(['Unknown ''Initialization'' option']);
end

% now check if tau is an array; if it is, it has to 
% have the same size as theta
if prod(size(tau)) > 1
   try,
      dummy = theta.*tau;
   catch,
      error(['Parameter tau has wrong dimensions; it should be scalar or size(theta)']),
   end
end
      
   

% if the true theta was given, check its size
if compute_mse & (size(true) ~= size(theta))  
   error(['Initial theta has incompatible size']); 
end

% initialize u, and v
u =  theta.*(theta >= 0);
v = -theta.*(theta <  0);

% define the indicator vector or matrix of nonzeros in theta
nz_theta = (theta ~= 0.0);

% Compute and store initial value of the objective function
resid =  y - A(theta);
prev_f = 0.5*(resid(:)'*resid(:)) + ...
         sum(tau(:).*u(:)) + sum(tau(:).*v(:));

times(1) = cputime - t0;
objective(1) = prev_f;

if compute_mse
   mses(1) = sum(sum((theta-true).^2));
end

% Compute the initial gradient and the useful 
% quantity resid_base
resid_base = y - resid;
temp = AT(resid_base);
term = temp - Aty;
gradu =   term + tau;
gradv =  -term + tau;

% projection and computation of search direction vector
du = max(u - gradu, 0.0) - u;
dv = max(v - gradv, 0.0) - v;

% initialization of alpha
alpha = 1/max(max(abs(du(:))),max(abs(dv(:))));
alphas(1) = alpha;

% control variable for the outer loop and iteration counter
% cont_outer = (norm(projected_gradient) > 1.e-5) ;

cont_outer = 1;
iter = 0;

while cont_outer

  theta_previous = theta;
  
  % compute gradient
  temp = AT(resid_base);
  
  term  =  temp - Aty;
  gradu =  term + tau;
  gradv = -term + tau;
  
  % projection and computation of search direction vector
  du = max(u - alpha*gradu, 0.0) - u;
  dv = max(v - alpha*gradv, 0.0) - v;
  dtheta = du-dv;
  old_u = u;
  old_v = v;
  
  % calculate useful matrix-vector product involving dtheta
  auv = A(dtheta);
  dGd = auv(:)'*auv(:);
  if (dGd > 0)
    lambda0 = - (gradu(:)'*du(:) + gradv(:)'*dv(:))/dGd;
    if (lambda0 <= 0)
      % if we get to here, there's something wrong
      alpha = alpha / 2; 
      fprintf(1,'alpha = %g\n',alpha);
      lambda = 0;
    else
     lambda = min(lambda0,1);
     u = old_u + lambda * du;
     v = old_v + lambda * dv;
     uvmin = min(u,v);
     u = u - uvmin; 
     v = v - uvmin; 
     theta = u - v;
    
    
    % update the nonzero indicator vector
     nz_theta_prev = nz_theta;
     nz_theta = (theta~=0.0);
     num_nz_theta = sum(nz_theta(:));
     num_changes_active = sum(nz_theta(:)~=nz_theta_prev(:));

     resid = y - resid_base - lambda*auv;
     f = 0.5*(resid(:)'*resid(:)) +  ...
        sum(tau(:).*u(:)) + sum(tau(:).*v(:));
     dd  = du(:)'*du(:) + dv(:)'*dv(:); 
     
     alpha = min(alphamax,max(alphamin,dd/dGd));
     resid_base = resid_base + lambda*auv;
    end  
    
  else %  (dGd < 0)
      
      % something wrong if we get to here
      fprintf(1,' dGd=%12.4e, something is wrong\n', dGd);
      alpha = alphamax;     
%       resid_base = resid_base + auv;
  end
  
  
% evaluate projected gradient quantity of Dai-Fletcher
%   projected_gradient = [min(gradu,0.0).*(u == 0.0) + gradu.*(u ~= 0.0);
%     min(gradv,0.0).*(v==0.0) + gradv.*(v~=0.0)];
  
  if ~stop_objective
     fprintf(1,'Iter=%4d, obj=%10.6e, alpha=%6.2e, lamda=%5.3f, nonzeros=%7d, chg=%6d\n',...
             iter,f,alpha,lambda,num_nz_theta,num_changes_active);
  else
     crit = abs(f-prev_f)/(prev_f*tolA);
     fprintf(1,'Iter=%4d, obj=%10.6e, alpha=%6.2e, lamda=%5.3f, stop crit.=%10.5e\n',...
             iter,f,alpha,lambda,crit);
  end
      
  
  iter = iter + 1;
  prev_f = f;
  objective(iter) = f;
  times(iter) = cputime-t0;
  alphas(iter) = alpha;
  
  if compute_mse
     err = true - theta;
     mses(iter) = (err(:)'*err(:));
  end
    
  % take no less than miniter and no more than maxiter iterations
  if iter<=miniter
     cont_outer = 1;
  else
    if ~stop_objective
       cont_outer = ((iter <= maxiter) & ...
                    ((num_changes_active > tolA*num_nz_theta)|...
                    (lambda<1.0)));
   else 
       cont_outer = ((iter <= maxiter) & ...
                    ((crit > 1)|(lambda<1.0)));
    end
           
  end
  
end % end of the main loop of the GPBB algorithm

% Printout results
fprintf(1,'\nFinished the main algorithm!\nResults:\n')
fprintf(1,'||A theta - y ||_2 = %10.3e\n',resid(:)'*resid(:))
fprintf(1,'||theta||_1 = %10.3e\n',sum(abs(theta(:))))
fprintf(1,'Objective function = %10.3e\n',f);
fprintf(1,'Number of non-zero components = %d\n\n',num_nz_theta)
% glb_nz(blkidx) = num_nz_theta;
% If the 'Debias' option is set to 1, we try to
% remove the bias from the l1 penalty, by applying CG to the 
% least-squares problem obtained by omitting the l1 term 
% and fixing the zero coefficients at zero.
if debias
    fprintf(1,'\n')
    fprintf(1,'Starting the debiasing phase...\n\n')

    theta_debias = theta;
    zeroind = (theta_debias~=0); 
    cont_debias_cg = 1;
    debias_start = iter;
  
    % calculate initial residual
    resid = A(theta_debias);
    resid = resid-y;
    resid_prev = eps*ones(size(resid));
    
    rvec = AT(resid);
  
    % mask out the zeros
    rvec = rvec .* zeroind;
    rTr_cg = rvec(:)'*rvec(:);
  
    % set convergence threshold for the residual || RW theta_debias - y ||_2
    tol_debias = tolD * (rvec(:)'*rvec(:));
  
    % initialize pvec
    pvec = -rvec;
  
    % main loop
    while cont_debias_cg
    
           % calculate A*p 
           RWpvec = A(pvec);      
           Apvec = AT(RWpvec);

           % mask out the zero terms
           Apvec = Apvec .* zeroind;
    
           % calculate alpha for CG
           if ((pvec(:)'* Apvec(:)) ~= 0)% Thong Do added this line to prevent 'divided by 0' error
               alpha_cg = rTr_cg / (pvec(:)'* Apvec(:));

               % take the step
               theta_debias = theta_debias + alpha_cg * pvec;
               resid = resid + alpha_cg * RWpvec;
               rvec  = rvec  + alpha_cg * Apvec;

               rTr_cg_plus = rvec(:)'*rvec(:);
               beta_cg = rTr_cg_plus / rTr_cg;
               pvec = -rvec + beta_cg * pvec;

               rTr_cg = rTr_cg_plus;
           end
    
           iter = iter+1;
    
           objective(iter) = 0.5*(resid(:)'*resid(:)) + ...
                             sum(tau(:).*abs(theta_debias(:)));
           times(iter) = cputime - t0;
    
           if compute_mse
              err = true - theta_debias;
              mses(iter) = (err(:)'*err(:));
           end

           if ~stop_objective
                fprintf(1,' Iter = %5d, debias resid = %13.8e, convergence = %8.3e\n', ...
     	                iter, resid(:)'*resid(:), rTr_cg / tol_debias);
               cont_debias_cg = ...
                   (iter-debias_start <= miniter_debias )| ...
                   ((rTr_cg > tol_debias) & ...
                   (iter-debias_start <= maxiter_debias));
           else
               err_prev = resid_prev(:)'*resid_prev(:);
               criteri = abs(resid(:)'*resid(:) - err_prev)/err_prev;
               fprintf(1,' Iter = %5d, debias resid = %13.8e, convergence = %8.3e\n', ...
	                iter, resid(:)'*resid(:), criteri / tol_debias);
               cont_debias_cg = ...
                   (iter-debias_start <= miniter_debias )| ...
                   ((criteri > tol_debias)  & ...
                   (iter-debias_start <= maxiter_debias));
           end
           resid_prev = resid;
    
    end
    fprintf(1,'\nFinished the debiasing phase!\nResults:\n')
    fprintf(1,'||A theta - y ||_2 = %10.3e\n',resid(:)'*resid(:))
    fprintf(1,'||theta||_1 = %10.3e\n',sum(abs(theta(:))))
    fprintf(1,'Objective function = %10.3e\n',f);
    nz = (theta_debias~=0.0);
    fprintf(1,'Number of non-zero components = %d\n\n',sum(nz(:)))
end

if compute_mse
   mses = mses/length(true(:));
end




