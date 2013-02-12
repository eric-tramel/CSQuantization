function [mse, xhat, x] = RunSampleRbp
% Function runs Relaxed Belief Propagation (RBP) reconstruction.
%
% RBP is an iterative estimation algorithm for estimation from linear
% measurements. In this script we demonstrate how it can be used to do
% estimation for compressive sensing. We consider the estimation problem
% where CS measurements are corrupt by AWGN and then quantized. We assume
% input vector X to be distributes i.i.d. from Gauss-Bernoulli distribution
% with sparsity ration rho. See paper for more details.
%
% "Optimal Quantization for Compressive Sensing under Message Passing
% Reconstruction" by U. Kamilov, V. Goyal, and S. Rangan <a href =
% "http://arxiv.org/abs/1102.4652">[arXiv]</a>.
%
% Ulugbek Kamilov, 2010, STIR @ MIT.
%
% See also CREATEQUANTIZER, QUANTIZEMEASUREMENTS, RECONSTRUCTRBP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define quantizer and create quantizer object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-binned quantizer levels (different from number of bins)
nLevels = 128; % Use even numbers for binned quantizers

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
% Estimation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length of the signal
n = 8000;

% Sparsity ratio
rho = 0.1;

% Measurement ratio (beta = n/m)
beta = 2;

% Number of measurements
m = ceil(n/beta);

% Noise variance
v = 0;

% Number of iteration of the algorithm
T = 100;

% Effective noise threshold
tol = 1e-3;

% Print messages to command-line during estimation
verbose = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the signal (Gaussian with proba. rho)
x = binornd(1, rho, n, 1) .* randn(n, 1) ./ sqrt(rho);

% Generate the measurement matrix
Phi = (1/sqrt(m)) .* randn(m, n);
Phi = Phi ./ repmat(sqrt(sum(Phi.^2)), m, 1);

% Obtain measurements
z = Phi*x;

% noise std. deviation
sigma = sqrt(v);

% Noise
eta = sigma .* randn(m, 1);

% Noisy measurements
s = z + eta;

% Quantize the measurements
y = quantizeMeasurements(s, q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run relaxed belief propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xhat, mse] = reconstructRbp(Phi, y, v, rho, q, x', T, tol, verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure
figure('Color', 'w');

% Plot quantized measurements
subplot(3, 1, 1);
plotQuantizedMeasurements(z, s, y, q);

% Plot MSE
subplot(3, 1, 2);
plot(1:length(mse), mse, 'LineWidth', 3);
set(gca, 'FontSize', 14);
grid on;
title('MSE of reconstruction');
xlabel('Iteration Number');
ylabel('MSE [dB]');
xlim([1, length(mse)]);

% Plot reconstruction
subplot(3, 1, 3);
stem(x, 'bo');
hold on;
stem(xhat, 'r.');
set(gca, 'FontSize', 14);
title('Reconstruction');
legend('signal', 'estimate');
xlim([1, n]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotQuantizedMeasurements(z, u, y, q)
% Helper function used to plot results of RunSampleRbp.

plot(z, 'm.');
hold on;
plot(u, 'r+');

m = length(z);

% Quantization decision boundaries
maxPlotEdge = q.x(end-1) + 5;

% Go through each quantized measurement
for i = 1:m
    % Obtain boundaries for this measurement
    boundaries = q.inverse{y(i)};
    
    % Number of boundaries
    nBoundaries = length(boundaries);
    
    % Number of cells
    nCells = nBoundaries/2;
    
    % Go through each cell and draw a rectangle
    for j = 1:nCells
        % Cell boundaries
        low = boundaries(2*j-1);
        high = boundaries(2*j);
        
        % Box height width and positions
        xbox = i-0.5;
        ybox = max(-maxPlotEdge, low);
        wbox = 1;
        hbox = min(maxPlotEdge, high) - ybox;
        
        % Draw rectangle
        rectangle('Position', [xbox, ybox, wbox, hbox], 'EdgeColor', 'c');
    end    
end

legend('Noiseless Measurement', 'Noisy Measurement');
set(gca, 'FontSize', 14);
title('Quantization of the measurements');
xlim([1, m]);