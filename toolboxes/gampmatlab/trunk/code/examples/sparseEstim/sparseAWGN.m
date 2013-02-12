% sparseAWGN:  Example of estimating a sparse vector with Gaussian noise.
%
% In this problem, x is a Bernoulli-Gaussian random vector, that we want to
% estimate from measurements of the form
%
%   y = A*x + w,
%
% where A is a random matrix and w is Gaussian noise.  This is a classical
% compressed sensing problem of estimating a sparse vector for random
% linear measurements.

% Set path to the main directory
addpath('../../main/');

% Parameters
nx = 200;         % Number of input components (dimension of x)
nz = 100;         % Number of output components (dimension of y)
sparseRat = 0.1;    % fraction of components of x that are non-zero
snr = 20;           % SNR in dB.

% Create a random sparse vector
xmean0 = 0;
xvar0 = 1;
x0 = normrnd(xmean0, sqrt(xvar0),nx,1); % a dense Gaussian vector
x = x0.*(rand(nx,1) < sparseRat);       % insert zeros

% Create a random measurement matrix
A = (1/sqrt(nx))*randn(nz,nx);

% Compute the noise level based on the specified SNR. Since the components
% of A have power 1/nx, the E|y(i)|^2 = E|x(j)|^2 = sparseRat.  
wvar = 10^(-0.1*snr)*xvar0*sparseRat;
w = normrnd(0, sqrt(wvar), nz, 1);
y = A*x + w;

% Generate input estimation class
% First, create an estimator for a Gaussian random variable (with no
% sparsity)
inputEst0 = AwgnEstimIn(xmean0, xvar0);

% Then, create an input estimator from inputEst0 corresponding to a random
% variable x that is zero with probability 1-sparseRat and has the
% distribution of x in inputEst0 with probability sparseRat.
inputEst = SparseScaEstim( inputEst0, sparseRat );

% Output estimation class:  Use the 
outputEst = AwgnEstimOut(y, wvar);

% Run the GAMP algorithm
opt = GampOpt();
[xhat, xvar] = gampEst(inputEst, outputEst, A, opt);

% Plot the results
h = stem([x xhat]);
grid on;
set(gca,'FontSize',16);
set(h, 'LineWidth', 2);
legend('True', 'GAMP Estimate');

% Display the MSE
mseGAMP = 20*log10( norm(x-xhat)/norm(x));
fprintf(1,'GAMP: MSE = %5.1f dB\n', mseGAMP);

