function e = stateEvolution(mzinit, rho, beta, v, quantizer, T, tol, c, nSamples, verbose)
% STATEEVOLUTION Performs state evolution for the setting in the paper
% (see reference below).
%
% e = stateEvolution(mzinit, T, rho, beta, v, quantizer, tol, c, nSamples, verbose)
%
% Input:
% - mzinit: Initialization of the algorithm (mzinit = beta * mxinit)
% In our simulations we have Gauss-Bernoulli X with mxinit = 1, where
% mxinit is variance of the elements of the signal vec.
% - rho: sparsity ratio of X.
% - beta = n/m: undersampling of the signal.
% - v: AWGN variance
% - quantizer: structure representing quantizer. To creat quantizers
% consult the function CREATEQUANTIZER.
% - T: maximum number of iterations of the algorithm (default T = 100).
% - tol: stopping criterion. If absolute change in MSE is smaller than tol
% stop (default tol = 0.0001).
% - c: parameter determining range for approximating infinite integrals.
% It indicates number of standard deviations to take for capturing most of
% the energy of the gaussian random variable.
% (default c = sqrt(2)*erfcinv(1e-12))
% - nSamples: number of samples for sampling PDFs (default nSamples = 2001)
% - verbose: Set to 1 to print messages to command line (default = 0)
%
% Output:
% - e: array of per iteration MSE predictions
%
% See "Optimal Quantization for Compressive Sensing under Message Passing
% Reconstruction" by U. Kamilov, V. Goyal, and S. Rangan <a href =
% "http://arxiv.org/abs/1102.4652">[arXiv]</a>.
%
% Ulugbek Kamilov, 2010, STIR @ MIT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10, verbose = 0; end;
if nargin < 9, nSamples = 2001; end;
if nargin < 8, c = sqrt(2)*erfcinv(1e-12); end;
if nargin < 7, tol = 0.0001; end;
if nargin < 6, T = 100; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize state evolution
mz = mzinit;

% Initialize output vector
e = zeros(T, 1);

% Iterate through the states
for t = 1:T
    % Perform recursion
    mq = eout(quantizer, mz, mzinit, v, c, nSamples);
    mz = beta * ein(mq, rho, c, nSamples);
    
    % Compute MSE at this iteration
    e(t) = 10*log10(mz/beta);
    
    % Print messages
    if(verbose)
        fprintf('mse(t = %d) = %.4f', t, e(t));
        if(t > 1)
            fprintf(', change = %.4f\n', abs(e(t)-e(t-1)));
        else
            fprintf('\n');
        end
    end
    
    % Verify stopping criteria
    if(t > 1 && abs(e(t) - e(t-1)) < tol)
        e = e(1:t);
        return;
    end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function e = ein(mu, rho, c, nSamples)
% Compute variable node update.

% Define temporary variable alpha
alpha = 1/rho;

% Variance of the wide gaussian
vg = alpha + mu;

% Rounded largest Q to consider
maxq = ceil(c*sqrt(vg));

% Sampling grid
sgrid = linspace(-maxq, maxq, nSamples);

% Gaussian of variance (alpha + mu) (= wide gaussian)
phi1 = normpdf(sgrid, 0, sqrt(vg));
phi1 = phi1/sum(phi1);

% Gaussian of variance (mu)
phi2 = normpdf(sgrid, 0, sqrt(mu));
phi2 = phi2/sum(phi2);

% Nominator
N = (sgrid.^2) .* (phi1.^2);

% Denominator
D = rho * phi1 + (1-rho) * phi2;
D(D<=0) = eps;

% Evaluate integral
I =  (vg.^(-2)) * sum(N ./ D);

% Result
e = 1 - I;

function e = eout(q, mz, mzinit, v, c, nSamples)
% Compute measurement node update.

% Number of quantization levels
N = q.nBins;

% Expected value
eval = 0;

% Effective noise variance
veff = mz + v;

% Go through each level
for y = 1:N
    % Extract quantization interval for this y
    qinv = q.inverse{y};
    
    % Compute integral over zhat
    I = evaluateZhatIntegral(mz, mzinit, v, qinv, nSamples, c);
    
    % Update expected value
    eval = eval + I;
end

% Result
e = (veff^2)/(veff - eval);

function y = evaluateZhatIntegral(mz, mzinit, v, qinv, nSamples, c)
% Computes the expected variance of Z ~ N(zhat, mz+v) over the interval
% given by qinv. Zhat itself is distributed according to N(0, mzinit-mz).
% Integral is replaced by sum via discretization of the grid.
%
% Input:
% - mz: Current state of SE
% - mzinit: Initial state of SE
% - v: AWGN variance
% - qinv: quantization invervals given as an array: [a1, b1, a2, b2, ...]
% - nSamples: number of samples to use to discretize zhat
% - c: c = sqrt(2)*erfcinv(epsilon), determines the bounds of integration
%       determined by acceptable error epsilon.
%
% Output:
% - yr = E{Var{Z | Y = y}} where y is the quantized measurement
%
% Ulugbek Kamilov, STIR, MIT, 2011.


% Standard deviation of the zhat
sigmaz = sqrt(mzinit - mz);

% Effective noise variance
veff = mz + v;

% If integral over delta function
if(sigmaz <= eps(0))
    [m0, m1, m2] = evaluateTotalMoment(0, veff, qinv);
    m0(m0 <= 0) = eps;
    y = m2 - (m1^2)/m0;
    return;
end

% Bounds of integration
maxz = ceil(c*sigmaz);

% Discrete grid
sgrid = linspace(-maxz, maxz, nSamples)';

% Normal pdf
gauspmf = normpdf(sgrid, 0, sigmaz);
gauspmf = gauspmf/sum(gauspmf);

% Compute moments
[m0, m1, m2] = evaluateTotalMoment(sgrid, veff, qinv);

% Make sure that m0 is not 0
m0(m0 <= 0) = eps;

% Compute function to integrate
f = m2 - (m1.^2)./m0;

% Integration
y = sum(f .* gauspmf);

function [m0, m1, m2] = evaluateTotalMoment(m, v, qinv)

% Extract bounds
av = qinv(1:2:end);
bv = qinv(2:2:end);

% Number of means
N = length(m);

% Number of intevals
M = length(av);

% Reshape arrays
av = repmat(av, N, 1);
bv = repmat(bv, N, 1);
m = repmat(m, 1, M);

% For each interval comput moment
[m0s, m1s, m2s] = evaluateMoment(m, v, av, bv);

% Retusn sum of the moments
m0 = sum(m0s, 2);
m1 = sum(m1s, 2);
m2 = sum(m2s, 2);

function [m0, m1, m2] = evaluateMoment(m, v, a, b)
% Evaluates moments 0,1,2 of a gaussian pdf in [a, b]
% a and b can be arrays of the same length

% Compute error functions
erfa = erf((a-m)/sqrt(2*v));
erfb = erf((b-m)/sqrt(2*v));

% Compute normal functions
gaua = normpdf(a, m, sqrt(v));
gaub = normpdf(b, m, sqrt(v));

% Evaluate moments
m0 = 0.5*(erfb - erfa);
m1= m .* m0 - v*(gaub - gaua);
m2= ((m.^2) + v).*m0 + v*((a+m).*gaua - (b+m).*gaub);