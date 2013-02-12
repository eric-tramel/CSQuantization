
% choose parameters
prior = 'ellp'; % e.g., 'bernoulli-gaussian' or 'laplacian' or 'bernoulli-ellp' ...
maxSum = 1;	% MAP (or max-sum) versus MMSE (or sum-product)
p1 = 0.2;	% "on" probability for bernoulli-* priors
rmax = 2;	% range of input mean
rvar = 0.25;	% input variance 

% build title
tit_str = prior;
if maxSum,
  tit_str = ['MAP ',tit_str];
else
  tit_str = ['MMSE ',tit_str];
end;
tit_str = [tit_str,', rvar=',num2str(rvar)];
if strfind(prior,'bernoulli'),
  tit_str = [tit_str,', p1=',num2str(p1)];
end;

% configure threshold
if strfind(prior,'gaussian')  
  xhat0 = 0;
  xvar0 = 1;
  inputEst0 = AwgnEstimIn(xhat0,xvar0,maxSum);
  tit_str = [tit_str,', xhat0=',num2str(xhat0),', xvar0=',num2str(xvar0)];

elseif strfind(prior,'laplacian')  
  if maxSum~=1, error('MMSE Laplacian not yet implemented!'); end;
  lambda = 0.5;
  inputEst0 = SoftThreshEstimIn(lambda);
  tit_str = [tit_str,', lambda=',num2str(lambda)];

elseif strfind(prior,'ellp')  
  if maxSum~=1, error('MMSE ell-p not yet implemented!'); end;
  lambda = 0.5;
  p = 1.0;
  inputEst0 = EllpEstimIn(lambda,p);
  tit_str = [tit_str,', p=',num2str(p),', lambda=',num2str(lambda)];

else
  error('Prior not recognized!')
end

% sparsify
if strfind(prior,'bernoulli'),
  inputEst = SparseScaEstim(inputEst0,p1,0,maxSum);
else 
  inputEst = inputEst0;
end

% compute input-output relationship
rhat = linspace(-rmax,rmax,1e5);
[xhat,xvar,val] = inputEst.estim(rhat,rvar*ones(size(rhat)));


% plot input-output relationship
clf;
subplot(121)
 handy = plot(rhat,xhat);
 set(handy,'Linewidth',2)
 hold on; plot(rmax*[-1,1],rmax*[-1,1],'r--'); hold off;
 xlabel('rhat');
 ylabel('xhat');
 axis('equal')
 grid on;
subplot(122)
 handy = plot(rhat,xvar);
 set(handy,'Linewidth',2)
 hold on; plot(rmax*[-1,1],rvar*[1,1],'r--'); hold off;
 xlabel('rhat');
 ylabel('xvar');
 grid on;
axes('Position',[0 0 1 1],'Visible','off');
 handy = text(0.5,0.96,tit_str);
 set(handy,'horizontalalignment','center')
