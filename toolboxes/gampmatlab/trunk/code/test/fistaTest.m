%Test code to test FISTA and compare to GAMP using soft thresholding
clear all
clc

% Set path
addpath('../main/');

%Handle random seed
defaultStream = RandStream.getDefaultStream;
if 0
    savedState = defaultStream.State;
    save random_state.mat savedState;
else
    load random_state.mat
end
defaultStream.State = savedState;

%% Setup and global options
%Specify problem parameters
N = 500;
del = 0.4;  	% ratio m/n (default 0.4)
rho = 0.2; 	% ratio k/m of sparsity to number of measurements (default 0.2)
SNR = 60;	% 
Atype = 'corr';	% type of matrix A: 'iid', 'corr', 'dft'
corr_tau = 0.9;	% column correlation for 'corr' type matrices (0 <= tau < 1)
dft_os = 4;	% column oversampling factor for 'dft' type matrices (dft_os=1,2,3,...)
shuffle_cols = 0; % shuffle (correlated) columns of A?

%Set options for GAMP
GAMP_options = GampOpt; %initialize the options object
%GAMP_options.verbose = 1;	
GAMP_options.nit = ceil(400/(1-corr_tau));% make sure we don't terminate too early [10000]
GAMP_options.tol = 1e-8;	% make sure we don't terminate too early [1e-8]
%GAMP_options.stepTol = 1e-10;  % terminate when step becomes very small [1e-10]
%GAMP_options.pvarMin = 1e-4;	% avoid numerical problems [1e-10]
%GAMP_options.xvarMin = 1e-10;	% avoid numerical problems [1e-10]
%GAMP_options.adaptStep = 1;	% definitely adapt the stepsize [1]
%GAMP_options.stepMax = 1;	% cap the stepsize at 1
GAMP_options.step = 1-corr_tau;	% starting small becomes essential as tau -> 1
GAMP_options.stepIncr = 1.1;	% important to increase stepsize when possible [1.1]
%GAMP_options.stepDecr = 0.4;	% decrease stepsize when needed [0.4]
%GAMP_options.stepWindow = 20;	% allow alg to take a few temporarily bad steps [20]
GAMP_options.pvarStep = 1;	% seems to work best with this on [1]
GAMP_options.varNorm = 0;	% seems to work best with this off [0] 
GAMP_options.uniformVariance = 1;% seems to work either way, so turn on for speed [1]

%Set FISTA options
FISTA_options = FistaOpt();
FISTA_options.tol = GAMP_options.tol; %same tolerance
FISTA_options.nit = 10000; %but often more iterations
FISTA_options.verbose = false;

% target value of Lambda/sqrt(muw) (set to <0 to use value corresponding to DMM output)
Lam_target = 1;		% normalized to be approximately independent of SNR

% DMM soft threshold parameter
alpha = 2.5; 			% 2.5 seems to be a decent choice

% ell-p norm parameter
p = 1.0;			% 1 for Lasso/BPDN


%% Generate the forward operator
%Derive sizes
M = ceil(del*N);
K = floor(rho*M);

%Avoid zero K
if K == 0
    K = 1;
end

%Do the AR model manually- actually faster when M is large
if strcmp(Atype,'iid'),
  A = randn(M,N) + 1j*randn(M,N);

elseif strcmp(Atype,'corr'),
  var0 = (1 - corr_tau)^2 / (1 - corr_tau^2); %initial variance
  A = zeros(M,N);
  A(:,1) = sqrt(var0/2)*(randn(M,1) + 1j*randn(M,1));
  for kk = 2:N
    A(:,kk) = corr_tau*A(:,kk-1) +  (1 - corr_tau)*sqrt(1/2)*(randn(M,1) + 1j*randn(M,1));
  end

elseif strcmp(Atype,'dft'),
  rows = randperm(N);
  A = fliplr(vander(exp(j*2*pi/N/dft_os*[0:N-1]')))';
  A = A(rows(1:M),:)/sqrt(M);

else
  error('unknown Atype')
end;

% scale columsn to unit-norm 
A = A*diag(1 ./ sqrt(diag(A'*A)));	

% shuffle columns due to avoid interacting with block-sparse signal?
if shuffle_cols, A = A(:,randperm(N)); end;


%% Generate the true signal -----------------------------------------------------------

%Determine true bits
truebits = false(N,1);
truebits(1:K) = true; %which bits are on

%Generate the signal
x = (sqrt(1/2)*(randn(N,1) + 1j*randn(N,1))) .* truebits;

%Generate the uncorrupted measurements
ztrue = A*x;

%% Output settings

%Output channel- we assume AWGN with mean thetahat and variance muw
muw = norm(ztrue)^2/M*10^(-SNR/10);

%Compute noisy output
y = ztrue + sqrt(muw/2)*(randn(M,1) + 1j*randn(M,1));


%% GAMP with Donoho/Maleki/Montanari-style thresholding ----------------

%Input channel
%inputEst = SoftThreshDMMEstimIn(alpha);
inputEst = EllpDMMEstimIn(alpha,p);

%Output channel
outputEst = CAwgnEstimOut(y, muw,true);

tic
[resGAMP,~,~,~,~,~,~,~, estHistGAMP] = ...
    gampEst(inputEst, outputEst, A, GAMP_options);
tGAMP_DMM = toc;

%Compute effective Lambda/sqrt(muw) 
LamDMM = 2*max(abs(A'*(A*(resGAMP) - y)))/sqrt(muw);
if Lam_target<0,
  Lam_target = LamDMM;
end;

%% MAP GAMP for Laplacian signal and Gaussian noise -------------------

%Input channel 
%inputEst2 = SoftThreshEstimIn(Lam_target);
inputEst2 = EllpEstimIn(Lam_target/2/sqrt(muw),p);

tic
GAMP_options.xvar0 = sqrt(muw)/Lam_target*norm(y)/norm(A,'fro');	
[resGAMP2,~,~,~,~,~,~,~, estHistGAMP2] = ...
    gampEst(inputEst2, outputEst, A, GAMP_options);
tGAMP_MAP = toc;
LamMAP = 2*max(abs(A'*(A*(resGAMP2) - y)))/sqrt(muw);
if (p==1)&(abs(LamMAP-Lam_target)/Lam_target>0.05),
  warning(['Lam_target=',num2str(Lam_target),' whereas Lam_MAP=',num2str(LamMAP)])
end;


%% FISTA ---------------------------------------------------------------


%Run FISTA
if (p~=1)&&(exist('FISTA_lp')==2)
  FISTAlp_options.max_iter = FISTA_options.nit;
  FISTAlp_options.thresh = FISTA_options.tol;
  tic
  [resFISTA] = FISTA_lp(...
  	zeros(size(x)),y,Lam_target*sqrt(muw),p,FISTAlp_options,A);
  tFISTA = toc;
  estHistFISTA.xhat = resFISTA*[1,1];
else,
  FISTA_options.lam = Lam_target*sqrt(muw); %set to match GAMP
  tic
  [resFISTA,estHistFISTA] = fistaEst(...
      zeros(size(x)),y,A,FISTA_options);
  tFISTA = toc;
  LamFISTA = 2*max(abs(A'*(A*(resFISTA) - y)))/sqrt(muw);
  if (p==1)&(abs(LamFISTA-Lam_target)/Lam_target>0.05),
    warning(['Lam_target=',num2str(Lam_target),' whereas Lam_FISTA=',num2str(LamFISTA)])
  end;
end



%% Show results

%Compute error values
errGAMP = zeros(size(estHistGAMP.xhat,2),1);
for kk = 1:length(errGAMP)
    errGAMP(kk) = norm(x - estHistGAMP.xhat(:,kk)) / norm(x);
end

errGAMP2 = zeros(size(estHistGAMP2.xhat,2),1);
for kk = 1:length(errGAMP2)
    errGAMP2(kk) = norm(x - estHistGAMP2.xhat(:,kk)) / norm(x);
end

errFISTA = zeros(size(estHistFISTA.xhat,2),1);
for kk = 1:length(errFISTA)
  errFISTA(kk) = norm(x - estHistFISTA.xhat(:,kk)) / norm(x);
end


%Show the results
figure(1); clf
subplot(211)
 plot(real(x),'ko')
 hold on
 plot(real(resFISTA),'r+')
 plot(real(resGAMP),'bx')
 plot(real(resGAMP2),'m^')
 hold off
 ylabel('real')
 legend('Truth','FISTA','GAMP-DMM','GAMP-MAP')
 grid on;
 if strcmp(Atype,'iid')
   title('iid Gaussian matrix')
 elseif strcmp(Atype,'corr')
   title(['Gaussian matrix with column correlation ',num2str(corr_tau)])
 elseif strcmp(Atype,'dft')
   title(['DFT matrix with oversampling factor ',num2str(dft_os)])
 end
subplot(212)
 plot(imag(x),'ko')
 hold on
 plot(imag(resFISTA),'r+')
 plot(imag(resGAMP),'bx')
 plot(imag(resGAMP2),'m^')
 hold off
 ylabel('imag')
 legend('Truth','FISTA','GAMP-DMM','GAMP-MAP')
 grid on;
 %axis([0 N -.2 5])

%Show convergence history
figure(2); clf
subplot(311)
 plot(20*log10(abs(errFISTA)),'r--')
 hold on
 plot(20*log10(abs(errGAMP)),'b-')
 plot(20*log10(abs(errGAMP2)),'m-')
 xlabel('iteration')
 ylabel('NMSE (dB)')
 legend('FISTA','GAMP-DMM','GAMP-MAP')
 title(['NMSE: '...
	'GAMP-DMM=' num2str(20*log10(errGAMP(end))),' dB'...
	',  GAMP-MAP=' num2str(20*log10(errGAMP2(end))),' dB'...
    	',  FISTA: ' num2str(20*log10(errFISTA(end))),'dB'...
    	%' dB;  Normalized Difference: 'num2str(norm(resFISTA - resGAMP)/norm(resGAMP))
    ]);
 axis('tight')
 grid on;
subplot(312)
 plot(estHistGAMP.step,'b-')
 hold on;
 plot(estHistGAMP2.step,'m-')
 hold off;
 legend('GAMP-DMM','GAMP-MAP')
 ylabel('stepsize');
 grid on;
 axis('tight'); axe=axis; axis([axe(1:2),0,1.05])
 title(['k-true=' num2str(sum(x~=0)),...
 	',  k-GAMP-DMM=' num2str(sum(resGAMP~=0)),...
	',  k-GAMP-MAP=' num2str(sum(resGAMP2~=0)),...
    	',  k-FISTA= ' num2str(sum(resFISTA~=0)),...
    ]);
subplot(313)
 semilogy(estHistGAMP.val,'b-')
 hold on;
 semilogy(estHistGAMP2.val,'m-')
 hold off;
 legend('GAMP-DMM','GAMP-MAP','Location','SouthEast')
 ylabel('val');
 grid on;
 axis([axe(1:2),min([estHistGAMP.val;estHistGAMP2.val]),max([estHistGAMP.val;estHistGAMP2.val])])
 if p==1,
   title(['lambda/sqrt(muw): '...
	'GAMP-DMM=' num2str(LamDMM,4),...
	',  target=' num2str(Lam_target,4),...
	',  GAMP-MAP=' num2str(LamMAP,4),...
    	',  FISTA= ' num2str(LamFISTA,4),...
    ]);
 else
   title(['p=',num2str(p)]);
 end;
