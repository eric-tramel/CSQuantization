function [mseSe, mseAmp] = RunSampleMonteCarlo
% Function performs Monte Carlo simulation to compare state evolution
% prediction to empirical performance of Approximate Message Passing (AMP)
% algorithm.
%
% AMP is an iterative estimation algorithm for estimation from linear
% measurements. In this function we demonstrate how it can be used to do
% estimation for compressive sensing. We consider the estimation problem
% where CS measurements are corrupt by AWGN and then quantized. We assume
% input vector X to be distributes i.i.d. from Gauss-Bernoulli distribution
% with sparsity ration rho. See paper for more details.
%
% NOTE: This script uses MATLAB's Parallel Computing Toolbox. For Monte
% Carlo simulations it uses default of 4 threads (c.f. matlabpool open 4;)
% Feel free to increase or reduce this number.
%
% "Optimal Quantization for Compressive Sensing under Message Passing
% Reconstruction" by U. Kamilov, V. Goyal, and S. Rangan <a href =
% "http://arxiv.org/abs/1102.4652">[arXiv]</a>.
%
% Ulugbek Kamilov, 2011, BIG @ EPFL.
%
% See also CREATEQUANTIZER, QUANTIZEMEASUREMENTS, RECONSTRUCTAMP, STATEEVOLUTION.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define quantizer and create quantizer object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-binned quantizer levels (different from number of bins)
nLevels = 32; % Use even numbers for binned quantizers

% This point replaces Inf
maxEdge = 1e2;

% Generate boundaries
x = [-maxEdge, linspace(-5, 5, nLevels-1), maxEdge];

% Binning rule - Regular quantizers
binTable = [(1:nLevels)', (1:nLevels)'];

% Binning rule - L=2 binned quantizers
%binTable = [(1:nLevels)', repmat((1:nLevels/2)', 2, 1)];

% Generate quantizer
q = createQuantizer(x, binTable);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sparsity ratio
rho = 0.1;

% Measurement ratio (beta = n/m)
beta = 2;

% Noise variance
v = 0.01;

% Number of iteration of the algorithm
T = 15;

% Stopping criterion. Stop if absolute change in MSE < tol. For this
% simulation set to 0.
tol = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Evolution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number samples of pdf (odd number)
nSamples = 2001;

% Small value for determining accuracy of numerical integrals
epsilon = 1e-14;

% Number of std. deviations to integrate
c = sqrt(2) * erfcinv(epsilon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run state evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mseSe = stateEvolution(beta, rho, beta, v, q, T, tol, c, nSamples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of monte carlo trials
C = 1000;

% Length of the signal
n = 2000;

% Number of measurements
m = ceil(n/beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run monte carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mseAmp = monteCarloAmp(C, T, n, rho, m, v, q, tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color', 'w');
subplot(1, 2, 1);
plot(1:T, mseSe, 'b-', 'LineWidth', 3);
hold on;
plot(1:T, median(mseAmp), 'r--', 'LineWidth', 3);
set(gca, 'FontSize', 14);
grid on;
title('State Evolution MSE prediction');
xlabel('Iteration Number');
ylabel('MSE [dB]');
xlim([1, T]);
legend('SE (pred)', sprintf('sim n = %d', n));

subplot(1, 2, 2);
min_mse = min(mseAmp(:, end));
max_mse = max(mseAmp(:, end));
nBins = 50;
bins = linspace(min_mse, max_mse, nBins);
chist = cumsum(hist(mseAmp(:, end), bins))/C;
plot([mseSe(end), mseSe(end)], [0, 1], 'b-', 'LineWidth', 3);
hold on;
plot(bins, chist, 'r--', 'LineWidth', 3);
set(gca, 'FontSize', 14);
grid on;
title('State Evolution MSE prediction');
xlabel('MSE [dB]');
ylabel('Cumul. Occurence');
legend('SE (pred)', sprintf('sim n = %d', n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mseAmp = monteCarloAmp(C, T, n, rho, m, v, q, tol)
% Performs Monte Carlo simulation of AMP estimation.
% Uses parallel computing toolbox of MATLAB. To disable it change 'parfor'
% to 'for', and remove lines with 'matlabpool'.

% Variable to hold the mse values
mseAmp = zeros(C, T);

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4;
end

parfor c = 1:C
    % Print to command line
    fprintf('Simulation %4.d out of %d\n', c, C);
    
    % Generate the signal (Gaussian with proba. rho)
    x = binornd(1, rho, n, 1) .* randn(n, 1) ./ sqrt(rho);
    
    % Generate the measurement matrix
    A = (1/sqrt(m)) .* randn(m, n);
    A = A ./ repmat(sqrt(sum(A.^2)), m, 1);
    
    % Obtain measurements
    z = A*x;
    
    % noise std. deviation
    sigma = sqrt(v);
    
    % Generate noise
    w = sigma .* randn(m, 1);
    
    % Noisy measurements
    u = z + w;
    
    % Quantize the measurements
    y = quantizeMeasurements(u, q);
    
    % Estimate x using linearized belief propagation
    [~, mseAmp(c, :)] = reconstructAmp(A, y, v, rho, q, x', T, tol);
end

matlabpool close;

% Remove errors
iFinite = ~isfinite(mseAmp);
iFinite = sum(iFinite, 2) == 0;
mseAmp = mseAmp(iFinite, :);