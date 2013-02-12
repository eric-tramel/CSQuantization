function mse = RunSampleStateEvolution
% Function runs State Evolution (SE) equations.
%
% SE predicts asymptotic MSE performance of Relaxed Belief Propagation. 
% The channel model this script considers is that of AWGN followed by a 
% Quantizer. We consider signals distributed according to Gauss-Bernoulli
% distribution with some parameter rho. See paper text for details.
%
% See "Optimal Quantization for Compressive Sensing under Message Passing
% Reconstruction" by U. Kamilov, V. Goyal, and S. Rangan <a href =
% "http://arxiv.org/abs/1102.4652">[arXiv]</a>.
%
% Ulugbek Kamilov, 2010, STIR @ MIT.
%
% See also STATEEVOLUTION, CREATEQUANTIZER.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define quantizer and create quantizer object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-binned quantizer levels (different from number of bins)
nLevels = 16;

% This point replaces Inf
maxEdge = 1e2;

% Generate boundaries
x = [-maxEdge, linspace(-5, 5, nLevels-1), maxEdge];

% Binning rule - Regular quantizers
% binTable = [(1:nLevels)', (1:nLevels)'];

% Binning rule - L=2 binned quantizers
binTable = [(1:nLevels)', repmat((1:nLevels/2)', 2, 1)];

% Generate quantizer
q = createQuantizer(x, binTable);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number samples of pdf (odd number)
nSamples = 2001;

% Small value for determining accuracy of numerical integrals
epsilon = 1e-12;

% Number of std. deviations to integrate
c = sqrt(2) * erfcinv(epsilon);

% Number of iteration of the algorithm
T = 100;

% Sparsity ratio
rho = 0.1;

% Measurement ratio (beta = n/m)
beta = 2;

% AWGN variance
v = 0.01;

% Stopping criterion
tol = 1e-4;

% Print messages to command line
verbose = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run state evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mse = stateEvolution(beta, rho, beta, v, q, T, tol, c, nSamples, verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color', 'w');
plot(1:length(mse), mse, 'LineWidth', 3);
set(gca, 'FontSize', 14);
grid on;
title('State Evolution MSE prediction');
xlabel('Iteration Number');
ylabel('MSE [dB]');
xlim([1, length(mse)]);