function [xhat, mse] = reconstructAmp(A, y, v, rho, quantizer, x, T, tol, verbose)
% RECONSTRUCTAMP Performs generalized approximate message passing
% estimation for the setting considered in the paper (see reference below).
%
% [xhat, mse] = reconstructAmp(A, y, v, rho, quantizer, x, T, tol, verbose)
%
% Input:
% - A: measurement matrix (m x n)
% - y: quantized measurements (typically indices) (m x 1)
% - v: AWGN variance
% - rho: sparsity ratio of the signal (rho = k/n)
% - quantizer: structure representing quantizer. To create quantizers
% consult the function CREATEQUANTIZER.
% - x: original signal (used to compute MSE performance) (1 x n)
% - T: maximum number of iterations of the algorithm (default T = 100)
% - tol: stopping criterion. If absolute change in MSE is smaller than tol
% stop (default tol = 0.0001)
% - verbose: print messages to command line
%
% Output:
% - xhat: final signal estimate (1 x n)
% - mse: array of per iteration MSE performance (T x 1)
%
% See "Optimal Quantization for Compressive Sensing under Message Passing
% Reconstruction" by U. Kamilov, V. Goyal, and S. Rangan <a href =
% "http://arxiv.org/abs/1102.4652">[arXiv]</a>.
%
% Ulugbek Kamilov, 2011, BIG @ EPFL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9, verbose = 0; end;
if nargin < 8, tol = 0.0001; end;
if nargin < 7, T = 100; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform main estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of measurements and dimension of the signal
[m, n] = size(A);

% AWGN standard deviation
sigma = sqrt(v);

% Initialize the estimates
xhat = zeros(1, n);
mx = ones(1, n);

% Initialize u (output update)
u = zeros(m, 1);

% Initialize MSE
mse = zeros(T, 1);

% Perform estimation
for t=1:T
    [u, mu] = outputUpdate(A, y, xhat, mx, u, sigma, quantizer);
    [xhat, mx] = inputUpdate(A, u, mu, xhat, rho);
    
    % Compute the empirical MSE and print it
    mse(t) = 10*log10(sum((xhat - x).^2)/n);
    
    % Print messages
    if(verbose)
        fprintf('mse(t = %d) = %.4f', t, mse(t));
        if(t > 1)
            fprintf(', change = %.4f\n', abs(mse(t)-mse(t-1)));
        else
            fprintf('\n');
        end
    end
    
    % Verify stopping criteria
    if(t > 1 && abs(mse(t) - mse(t-1)) < tol)
        mse = mse(1:t);
        return;
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xhat, mx] = inputUpdate(A, u, mu, xhat, rho)
% Performs input update

% Signal dimension
n = size(A, 2);

% Combine inputs
mq = 1 ./ sum((A.^2) .* repmat(mu, 1, n), 1);
q = xhat + mq .* sum(A .* repmat(u, 1, n), 1);

% Perform estimation
[xhat, mx] =  gin(q, mq, rho);

function [fin, ein] = gin(qv, mqv, rho)
% Compute non-linear functions Fin and Ein for each variable node
% These functions represent E{X|Q=q} and Var{X|Q=q} respectively, where
% Q=X+V, with X ~ p(x) with p(x) Gauss-Bernoulli of parameter rho, and
% V~Norm with mean 0 and variance mqv.
%
% Note: for very low variances this function might have numerical
% problems. To solve this issue PDF evaluations should be done in log
% domain, as it is commonly done in machine learning [e.g. instead of
% evaluating p(x) evaluate log(p(x)) then do operations in log domain
% directly].

% Signal length
n = length(qv);

% Temporary constant
alpha = 1/rho;

% Initialize outputs
fin = zeros(1, n);
ein = zeros(1, n);

% Perform estimation
for i = 1:n
    
    % Measurement and variance of effective AWGN
    q = qv(i);
    mu = mqv(i);
    
    % Compute the likelyhood of the measurement Q=q
    phi1 = normpdf(q, 0, sqrt(alpha+mu));
    phi2 = normpdf(q, 0, sqrt(mu));
    pq = rho*phi1 + (1-rho)*phi2;
    
    % Compute E{X|Q=q}
    t1 = phi1/((alpha+mu)*pq);
    fin(i) = q * t1;
    
    % Compute Var{X|Q=q}
    t2 = (alpha*q*q)/(alpha+mu);
    ein(i) = t1 * (mu + t2 - q*q*t1);
end

function [u, mu] = outputUpdate(A, y, xhat, mx, u, sigma, quantizer)
% Performs output node update.
%
% NOTE: This function can potentially run into numerical erros. This is due
% to the sub-function evaluateTotalMoment, which performs integration 
% of a gaussian in some integral given by quantizer boundaries. In case
% when this inteval is far from the mean of the normal and the normal has a
% small variance moments might result in 0, although in reality they
% represent some small values, ratio of which is definetely non-zero.

% length of the signal to estimate
m = size(A, 1);

% Combine inputs
mz = sum((A.^2) .* repmat(mx, m, 1), 2);
z = sum(A .* repmat(xhat, m, 1), 2) - (mz .* u);

% Total effective noise (AWGN + estiamtion)
mtv = mz + (sigma^2);

% Initialize outputs
u = zeros(m, 1);
mu = zeros(m, 1);

% Temporary array for measurements
yt = y;

% For each measurement node
while(~isempty(yt))
    % Pick the first non-considered measurement
    mes = yt(1);
    
    % Indices of the same measurements in y
    ix = y == mes;
    
    % Quantization region
    qinv = quantizer.inverse{mes};
    
    % Mean and variance
    zhat = z(ix);
    v = mtv(ix);
    
    % Compute moments (Note: might result in 0 moments due to num. errs)
    [m0, m1, m2] = evaluateTotalMoment(zhat, v, qinv);
    m0(m0<=0) = eps;
    
    % Set D1 and D2
    u(ix) = (m1./m0 - zhat)./v;
    mu(ix) = (1 - ((m2./m0) - (m1./m0).^2)./v)./v;
    
    % Modify yt
    yt = yt(yt ~= mes);
end

function [m0, m1, m2] = evaluateTotalMoment(m, v, qinv)
% Evaluates moment for each subinterval of qinv and combines. This function
% vectorizes values to make the evaluation fast.

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
v = repmat(v, 1, M);

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
erfa = erf((a-m)./sqrt(2*v));
erfb = erf((b-m)./sqrt(2*v));

% Compute normal functions
gaua = normpdf(a, m, sqrt(v));
gaub = normpdf(b, m, sqrt(v));

% Evaluate moments
m0 = 0.5*(erfb - erfa);
m1= m .* m0 - v.*(gaub - gaua);
m2= ((m.^2) + v).*m0 + v.*((a+m).*gaua - (b+m).*gaub);