function results = RunQuantizerOptimization
% Function performs quantizer optimization for given set of rates. All
% parameters are defined within the function.
% Optimization relies on State Evolution (SE) equations.
%
% results = RunQuantizerOptimization();
%
% Output:
% - results: structure containing six cell arrays: 
%   (1) results.quniform - uniform quantizer.
%   (2) results.qlloyds - Lloyd algorithm optimized quantizer
%   (3) results.qrouniform - regular uniform quantizer optimized via SE.
%   (4) results.qroptimal - optimal regular quantizer obtained via SE.
%   (5) results.qbuniform - binned uniform quantizer optimized via SE.
%   (6) results.qboptimal - optimal binned quantizer obtained via SE.
%
% Elements of the cell array are the results of optimization for each rate
% Rz for compressive sesning measurements (bits/measurement). Each element
% is another structure with two fields: (i) struct.mse - reconstruction
% error for the best quantizer, (ii) struct.quantizer - quantizer structure
% containing the best quantizer for this rate. For more information on
% quantizer structure see CREATEQUANTIZER.
%
% Additionally saveResults parameter can be set to 1, to save the results
% on hard drive as a MAT file and log the process as TXT file. See code for
% details.
%
% This file is mostly self contained, i.e. besides standard MATLAB 
% functions it relies only on two external functions STATEEVOLUTION and 
% CREATEQUANTIZER.
%
% File can be organized into 4 parts grouped together:
%   (a) this function (RUNQUANTIZEROPTIMIZATION) which sets all up.
%   (b) Wrapper functions for optimization (e.g. OPTIMIZEUNIF, ...)
%   (c) Helper functions (e.g. WRAPPERSE, ...)
%   (d) Lloyd max code (DESIGNGAUSSUNIFORMQUANTIZER,...)
%
% The part (d) was developed by Peter Kabal and dowloaded online. It
% allowys free redistribution and modification of the code as long as the
% copyright notice is included. Copyright notice can be found in the header
% commends of part (d) as well as README file.
%
% See "Optimal Quantization for Compressive Sensing under Message Passing
% Reconstruction" by U. Kamilov, V. Goyal, and S. Rangan <a href =
% "http://arxiv.org/abs/1102.4652">[arXiv]</a>.
%
% Ulugbek Kamilov, 2010, STIR @ MIT.
%
% See also STATEEVOLUTION, CREATEQUANTIZER.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimization algorithm (Use either 'interior-point' or 'active-set')
algo = 'interior-point';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SE input Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of iterations of the SE algorithm
T = 100;

% Sparsity ratio
rho = 0.1;

% Number of samples to use for PDF discretization
nSamples = 2001;

% Small value used as round off parameter
epsilon = 1e-12;

% Number of std. dev.
c = sqrt(2) * erfcinv(epsilon);

% AWGN variance (can be set to zero)
v = 0;

% Stopping criteria for state evolution
tol = 0.001;

% Undersampling (beta = n/m)
beta = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantizer parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximal value for boundary points (replaces infinity)
maxEdge = 1e2;

% Rates of measurements to consider (bits/measurement)
Rz = [2, 3, 4, 5];

% Number of quantizer outputs per bin
Ls = [2, 3];

% Number of Rz
nRz = length(Rz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for saving the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save and log results to disk?
saveResults = 0;

% If yes create log file and mat file
if(saveResults)
    timestamp = datestr(now, 'yyyy_mmm_dd_HH_MM');
    fileName = sprintf('RunQuantizerOptimization_%s.mat', timestamp);
    textFile = sprintf('RunQuantizerOptimization_%s.txt', timestamp);
    
    % Open text file logging
    fid = fopen(textFile, 'w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for running in parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save and log results to disk?
runParallel = 1;

% Is another parallel job open
isOpen = matlabpool('size') > 0;

% Number of cores to use
nCores = 4;

% Prepare for parallel processing
if(runParallel && ~isOpen) % checking to see if my pool is already open
    matlabpool('open', nCores);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save outputs of optimization in a structure with
% struct.mse: MSE prediction of SE
% struct.quantizer: quantizer struture
%
% Quantizer structure is defined as follows
% quantizer.nLevels: number of levels of underlying reg. quantizer
% quantizer.nBins: number of bins
% quantizer.inverse: Q^{-1}(y) function
% quantizer.x: quantizer edges

% Uniform quantizer (regular)
quniform = cell(nRz, 1);

% Lloyds quantizer (regular)
qlloyds = cell(nRz, 1);

% Optimal regular uniform quantizer
qrouniform = cell(nRz, 1);

% Optimal regular quantizer
qroptimal = cell(nRz, 1);

% Optimal binned uniform quantizer
qbuniform = cell(nRz, 1);

% Optimal binned quantizer
qboptimal = cell(nRz, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each rate find optimal quantizers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Regular
for i = 1:nRz
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize uniform quantizer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quniform{i} = optimizeUnif(Rz(i), beta, T, rho, v, maxEdge, tol, c, nSamples);
    
    % Print to command line
    fprintf('[Rz = %.2f] uniform quantizer achieving MSE = %.4f\n', Rz(i), quniform{i}.mse);
    
    % Log and save the result
    if(saveResults)
        fprintf(fid, '[Rz = %.2f] uniform quantizer achieving MSE = %.4f\r\n', Rz(i), quniform{i}.mse);
        save(fileName);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize lloyds quantizer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qlloyds{i} = optimizeLloyd(Rz(i), beta, T, rho, v, maxEdge, tol, c, nSamples);
    
    % Print to command line
    fprintf('[Rz = %.2f] Lloyd quantizer achieving MSE = %.4f\n', Rz(i), qlloyds{i}.mse);
    
    % Log and save the result
    if(saveResults)
        fprintf(fid, '[Rz = %.2f] Lloyd quantizer achieving MSE = %.4f\r\n', Rz(i), qlloyds{i}.mse);
        save(fileName);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize regular optimal uniform quantizer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization
    x0 = 3;
    
    qrouniform{i} = optimizeRegOptUnif(x0, Rz(i), beta, T, rho, v, maxEdge, tol, c, nSamples, algo);
    
    % Print to command line
    fprintf('[Rz = %.2f] optimal regular uniform quantizer achieving MSE = %.4f\n', Rz(i), qrouniform{i}.mse);
    
    % Log and save the result
    if(saveResults)
        fprintf(fid, '[Rz = %.2f] optimal regular uniform quantizer achieving MSE = %.4f\r\n', Rz(i), qrouniform{i}.mse);
        save(fileName);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize regular optimal quantizer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qroptimal{i} = optimizeRegOptimal(qrouniform{i}.quantizer.x, Rz(i), beta, T, rho, v, maxEdge, tol, c, nSamples, algo);
    
    % Print to command line
    fprintf('[Rz = %.2f] optimal regular quantizer achieving MSE = %.4f\n', Rz(i), qroptimal{i}.mse);
    
    % Log and save the result
    if(saveResults)
        fprintf(fid, '[Rz = %.2f] optimal regular quantizer achieving MSE = %.4f\r\n', Rz(i), qroptimal{i}.mse);
        save(fileName);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize uniform binned quantizer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization
    x0 = 5;
    
    qbuniform{i} = optimizeBinOptUnif(x0, Rz(i), Ls, beta, T, rho, v, maxEdge, tol, c, nSamples, algo);
    
    % Print to command line
    fprintf('[Rz = %.2f] optimal binned uniform  quantizer achieving MSE = %.4f\n', Rz(i), qbuniform{i}.mse);
    
    % Log and save the result
    if(saveResults)
        fprintf(fid, '[Rz = %.2f] optimal binned uniform quantizer achieving MSE = %.4f\r\n', Rz(i), qbuniform{i}.mse);
        save(fileName);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize binned quantizer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization range
    x0 = [0.05, 33];
    
    qboptimal{i} = optimizeBinOptimal(x0, Rz(i), Ls, beta, T, rho, v, maxEdge, tol, c, nSamples, algo);
    
    % Print to command line
    fprintf('[Rz = %.2f] optimal binned quantizer achieving MSE = %.4f\n', Rz(i), qboptimal{i}.mse);
    
    % Log and save the result
    if(saveResults)
        fprintf(fid, '[Rz = %.2f] optimal binned quantizer achieving MSE = %.4f\r\n', Rz(i), qboptimal{i}.mse);
        save(fileName);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close the file
if(saveResults)
    fclose(fid);
end

% Close the parallel computing
if(runParallel)
    matlabpool('close');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine results into single structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.quniform = quniform;
results.qlloyds = qlloyds;
results.qrouniform = qrouniform;
results.qroptimal = qroptimal;
results.qbuniform = qbuniform;
results.qboptimal = qboptimal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Define MSE arrays
mse_quniform = zeros(nRz, 1);
mse_qlloyds = zeros(nRz, 1);
mse_qrouniform = zeros(nRz, 1);
mse_qroptimal = zeros(nRz, 1);
mse_qbuniform = zeros(nRz, 1);
mse_qboptimal = zeros(nRz, 1);

% Populate array
for i = 1:nRz
    mse_quniform(i) = quniform{i}.mse;
    mse_qlloyds(i) = qlloyds{i}.mse;
    mse_qrouniform(i) = qrouniform{i}.mse;
    mse_qroptimal(i) = qroptimal{i}.mse;
    mse_qbuniform(i) = qbuniform{i}.mse;
    mse_qboptimal(i) = qboptimal{i}.mse;
end

figure('Color', 'White');
plot(Rz, mse_qlloyds, '--', Rz, mse_qroptimal, '-.',Rz, mse_qboptimal, 'LineWidth', 3);
legend('Lloyds', 'Optimal', 'Binned');
grid on;
set(gca, 'FontSize', 14);
xlabel('rate (bits/measurement)');
ylabel('MSE [dB]');
title(sprintf('Meas. Rate vs. Distortion: beta = %.2f', beta));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrapper optimization functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q = optimizeUnif(Rz, beta, T, rho, v, maxEdge, tol, c, nSamples)

% Get number of levels
nLevels = 2^Rz;

% Binning table for regular
binTable = repmat((1:nLevels)', 1, 2);

% Generate function pointers (this is for Lloyd-Max algorithm)
fpdf = gaussianPdf(beta, v);

% Run lloyd max
[~, quantx] = designGaussUniformQuantizer (nLevels, fpdf, 1);

% Perform state evolution of Lloyd-Max partitions
[mse, quantizer] = wrapperSe(quantx, beta, T, rho, beta, v, 'standard', binTable, maxEdge, tol, c, nSamples);

% Create structure
q.mse = mse;
q.quantizer = quantizer;

function q = optimizeLloyd(Rz, beta, T, rho, v, maxEdge, tol, c, nSamples)

% Get number of levels
nLevels = 2^Rz;

% Binning table for regular
binTable = repmat((1:nLevels)', 1, 2);

% Generate function pointers (this is for Lloyd-Max algorithm)
fpdf = gaussianPdf(beta, v);

% Run lloyd max
[~, quantx] = designGaussLloydQuantizer(nLevels, fpdf, 1);

% Perform state evolution of Lloyd-Max partitions
[mse, quantizer] = wrapperSe(quantx, beta, T, rho, beta, v, 'standard', binTable, maxEdge, tol, c, nSamples);

% Create structure
q.mse = mse;
q.quantizer = quantizer;

function q = optimizeRegOptUnif(x0, Rz, beta, T, rho, v, maxEdge, tol, c, nSamples, algo)

% Get number of levels
nLevels = 2^Rz;

% Options for optimization with fmincon
options = optimset(...
    'Algorithm', algo,...
    'Diagnostics', 'on',...
    'Display', 'iter-detailed',...
    'TolCon', 0,...
    'UseParallel','always');

% Binning table for regular
binTable = repmat((1:nLevels)', 1, 2);

% Get optimal boundaries
x = fmincon(@(x) wrapperSe(...
    x, beta, T, rho, beta, v, 'u', binTable, maxEdge, tol, c, nSamples),...
    x0, [], [], [], [], 0, maxEdge, [], options);

% Evaluate State Evolution
[mse, quantizer] = wrapperSe(x, beta, T, rho, beta, v, 'u', binTable, maxEdge, tol, c, nSamples);

% Create structure
q.mse = mse;
q.quantizer = quantizer;

function q = optimizeRegOptimal(x0, Rz, beta, T, rho, v, maxEdge, tol, c, nSamples, algo)

% Get number of levels
nLevels = 2^Rz;

% Number of edges
nEdges = nLevels - 1;

% Number of variables
nVariables = floor(0.5*nEdges);

% Options for optimization with fmincon
options = optimset(...
    'Algorithm', algo,...
    'Diagnostics', 'on',...
    'Display', 'iter-detailed',...
    'TolCon', 0,...
    'UseParallel','always');

% Binning table for regular
binTable = repmat((1:nLevels)', 1, 2);

% Initialize to optimal regular uniform
x0 = x0(end-1);
x0 = linspace(0, x0, nVariables+1);
x0 = x0(2:end);

% Construct constraints (Regular quantization)
nConstraints = nVariables-1;
A = zeros(nConstraints, nVariables);
for j = 1:nConstraints
    A(j, :) = [zeros(1, j-1), 1, -1, zeros(1, nConstraints-j)];
end
b = zeros(nConstraints, 1);

% Get optimal boundaries
x = fmincon(@(x) wrapperSe(...
    x, beta, T, rho, beta, v, 's', binTable, maxEdge, tol, c, nSamples),...
    x0, A, b, [], [], 0, maxEdge, [], options);

% Evaluate State Evolution
[mse, quantizer] = wrapperSe(x, beta, T, rho, beta, v, 's', binTable, maxEdge, tol, c, nSamples);

% Create structure
q.mse = mse;
q.quantizer = quantizer;

function q = optimizeBinOptUnif(x0, Rz, Ls, beta, T, rho, v, maxEdge, tol, c, nSamples, algo)

% Number of bins to consider
nLs = length(Ls);

% Number of bins
B = 2^Rz;

optmse = 1e6;

% Go through each bin
for i = 1:nLs
    % Number of levels per bin
    L = Ls(i);
    
    % Quantizer levels (non-binned)
    nLevels = L * B;
    
    % Binning rule
    binTable = [(1:nLevels)', repmat((1:B)', L, 1)];
    
    % Number of variables
    nVariables = floor(0.5*(nLevels-1));
    
    % Options for optimization with patternsearch
    psoptions = psoptimset(...
        'CompleteSearch','off',...
        'Display', 'diagnose',...
        'UseParallel', 'never',...
        'MaxIter', 20*nVariables,...
        'TimeLimit', 3 * 60 * 60,...
        'TolCon', 0,...
        'TolFun', 1e-4,...
        'Vectorized', 'off');
    
    % Get optimal boundaries
    x0 = patternsearch(...
        @(x) wrapperSe(...
        x, beta, T, rho, beta, v, 'u', binTable, maxEdge, tol, c, nSamples),...
        x0, [], [], [], [], 0, maxEdge, [], psoptions);
    
    % Options for optimization with fmincon
    options = optimset(...
        'Algorithm', algo,...
        'Diagnostics', 'on',...
        'Display', 'iter-detailed',...
        'TolCon', 0,...
        'UseParallel','always');
    
    % Perform fmincon
    [x, mse] = fmincon(@(x) wrapperSe(...
        x, beta, T, rho, beta, v, 'u', binTable, maxEdge, tol, c, nSamples),...
        x0, [], [], [], [], 0, maxEdge, [], options);
    
    % Check optimality of MSE
    if(mse < optmse)
        optmse = mse;
        optedges = x;
        optBinTable = binTable;
    end
end

% Evaluate State Evolution
[mse, quantizer] = wrapperSe(optedges, beta, T, rho, beta, v, 'u', optBinTable, maxEdge, tol, c, nSamples);

% Create structure
q.mse = mse;
q.quantizer = quantizer;

function q = optimizeBinOptimal(x0, Rz, Ls, beta, T, rho, v, maxEdge, tol, c, nSamples, algo)

% Number of bins to consider
nLs = length(Ls);

% Number of bins
B = 2^Rz;

optmse = 1e6;

% Go through each bin
for i = 1:nLs
    % Number of levels per bin
    L = Ls(i);
    
    % Quantizer levels (non-binned)
    nLevels = L * B;
    
    % Binning rule
    binTable = [(1:nLevels)', repmat((1:B)', L, 1)];
    
    % Number of variables
    nVariables = floor(0.5*(nLevels-1));
    
    % Options for optimization with patternsearch
    psoptions = psoptimset(...
        'CompleteSearch','off',...
        'Display', 'diagnose',...
        'UseParallel', 'never',...
        'MaxIter', 5*nVariables,...
        'TimeLimit', 3 * 60 * 60,...
        'TolCon', 0,...
        'TolFun', 1e-4,...
        'Vectorized', 'off');
    
    % Initial points
    xm0 = linspace(x0(1), x0(2), nVariables);
    
    % Construct constraints (Regular quantization)
    nConstraints = nVariables-1;
    A = zeros(nConstraints, nVariables);
    for j = 1:nConstraints
        A(j, :) = [zeros(1, j-1), 1, -1, zeros(1, nConstraints-j)];
    end
    b = zeros(nConstraints, 1);
    
    % Perform pattern search
    xm0 = patternsearch(...
        @(x) wrapperSe(x, beta, T, rho, beta, v, 's', binTable, maxEdge, tol, c, nSamples),...
        xm0, A, b, [], [], 0, maxEdge, [], psoptions);
    
    % Depending on algo set number of function evaluations
    if(strcmp(algo, 'active-set'))
        nEvals = 200 * nVariables;
        nIter = 800;
    else
        nEvals = 6000;
        nIter = 2000;
    end
    
    % Options for optimization with fmincon
    options = optimset(...
        'Algorithm', algo,...
        'Diagnostics', 'on',...
        'Display', 'iter-detailed',...
        'TolCon', 0,...
        'MaxFunEvals', nEvals,...
        'MaxIter', nIter,...
        'UseParallel','always');
    
    % Construct constraints (Regular quantization)
    nConstraints = nVariables-1;
    A = zeros(nConstraints, nVariables);
    for j = 1:nConstraints
        A(j, :) = [zeros(1, j-1), 1, -1, zeros(1, nConstraints-j)];
    end
    b = zeros(nConstraints, 1);
    
    % Perform fmincon
    [x, mse] = fmincon(@(x) wrapperSe(...
        x, beta, T, rho, beta, v, 's', binTable, maxEdge, tol, c, nSamples),...
        xm0, A, b, [], [], 0, maxEdge, [], options);
    
    % Check optimality of MSE
    if(mse < optmse)
        optmse = mse;
        optedges = x;
        optBinTable = binTable;
    end
end


% Evaluate State Evolution
[mse, quantizer] = wrapperSe(optedges, beta, T, rho, beta, v, 's', optBinTable, maxEdge, tol, c, nSamples);

% Create structure
q.mse = mse;
q.quantizer = quantizer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e, q] = wrapperSe(quantx, mzinit, T, rho, beta, v, type, binTable, maxEdge, tol, c, nSamples)

% Number of levels (total)
nLevels = size(binTable, 1);

% Number of edges
nEdges = nLevels - 1;

% Type of quantizer to optimize
% 1: uniform - generate nEdges points between [-quantizer, quantizer].
% Additionally prefix and postfix with +- maxEdge.
% 2: symmetric - flip the quantizer to generate symmetric boundaries around
% 0. Additioanlly add +- maxEdge.
% 3: complete - just pad with +- maxEdge.
switch type
    case {'uniform', 'u'}
        quantx = [-maxEdge, linspace(-quantx, quantx, nEdges), maxEdge];
    case {'symmetric', 's'}
        if(mod(nEdges, 2) == 0)
            quantx = [-maxEdge, -fliplr(quantx), quantx, maxEdge];
        else
            quantx = [-maxEdge, -fliplr(quantx), 0, quantx, maxEdge];
        end
    otherwise
        quantx = [-maxEdge, quantx, maxEdge];
end

% Generate quantizer
q = createQuantizer(quantx, binTable);

% Check if the boundaries are valid
if(~all(sort(quantx) == quantx))
    e = 1e6;
    return;
end

% Return result of state evolution
e = stateEvolution(mzinit, rho, beta, v, q, T, tol, c, nSamples);
e = e(end);

function FPDF = gaussianPdf(beta, v)
% Helper function that returns handles to 3 functions used to compute the
% area, the mean, and the second moment of gaussian distribution over some
% region of real axis. The aussian has mean 0 and variance (beta+v)

FPDF = {@(a, b)(gausArea(a, b, beta, v)),...
    @(a, b)(gausMean(a, b, beta, v)),...
    @(a, b)(gausVar(a, b, beta, v))};

function f = gausArea(a, b, beta, v)
f = normcdf(b, 0, sqrt(beta+v)) - normcdf(a, 0, sqrt(beta+v));

function f = gausMean(a, b, beta, v)
f = quadgk(@(x)(x .* normpdf(x, 0, sqrt(beta+v))), a, b);

function f = gausVar(a, b, beta, v)
f = quadgk(@(x)((x.^2) .* normpdf(x, 0, sqrt(beta+v))), a, b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Traditional optimization routines (Lloyds and Uniform)
%
% Main functions: 
% - designGaussUniformQuantizer
% - designGaussLloydQuantizer
%
% Downloaded from
% http://www.mathworks.com/matlabcentral/fileexchange/24333.
%
% Uploader allows redistribution and modification under condition that the
% following copyright notice is included.
%
% Copyright (c) 2009, Peter Kabal
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the McGill University nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Yq, Xq, MSE, Entropy, SNRdB, sdV] = designGaussUniformQuantizer (nLevels, FPDF, Sym)
% Find a uniform minimum mean square error quantizer.
%
% This subroutine searches for a uniform quantizer which minimizes the
% mean square quantization error. It returns the quantizer decision levels
% and output levels. A uniform quantizer is by our definition, one which
% has equally spaced output levels. Given those output levels, the
% MMSE quantizer has decision levels which fall midway between the output
% levels, and so are also equally spaced.
% nLevels - Number of quantizer output levels
% FPDF - Cell array of function handles {Farea, Fmean, Fvar}
% Farea - Pointer to a function to calculate the integral of a probability
%   density function in an interval.
%                 b
%   Farea(a,b) = Int p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fmean - Pointer to a function to calculate the mean of a probability
%   density function in an interval.
%                 b
%   Fmean(a,b) = Int x p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fvar - Pointer to a function to calculate the second moment of a
%   probability density function in an interval.
%                b
%   Fvar(a,b) = Int x^2 p(x) dx,
%                a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Sym - Logical flag. If true, the quantizer is required to be symmetric
%   about the mean. For symmetric quantizers, the probability density is
%   assumed to be symmetric about its mean.
%
% Yq - nLevels quantizer output levels (in ascending order)
% Xq - nLevels-1 quantizer decision levels (in ascending order)
% MSE - Resulting mean square error
% Entropy - Resulting entropy (bits)
% SNRdB - Resulting signal-to-noise ratio dB

Farea = FPDF{1};
Fmean = FPDF{2};

% Parameter
TolP = 1e-8;    % Symmetry check tolerance

% Consistency check
if (Sym)
  Xmean = feval(Fmean, -Inf, Inf);
  PH = feval(Farea, Xmean, Inf);
  if (abs(PH - 0.5) > TolP)
    fprintf('QuantUnif: Warning, PDF is not symmetric\n');
  end
end

% Find a minimum mean square error uniform quantizer
Yq = QuantUOpt(nLevels, FPDF, Sym);

% Generate the quantizer decision levels
if (nLevels > 1)
  Xq = 0.5 * (Yq(1:end-1) + Yq(2:end)); 
else
  Xq = [];
end

% Calculate the MSE and entropy
if (nargout > 2)
  [MSE, Entropy, SNRdB, sdV] = QuantSNR(Yq, Xq, FPDF);
end

return

function [Yq, Xq, MSE, Entropy, SNRdB] = designGaussLloydQuantizer (nLevels, FPDF, Sym)
% Find a non-uniform minimum mean square error quantizer.
%
% This subroutine searches for a quantizer which minimizes the mean
% square quantization error. It returns the quantizer decision levels
% and output levels that satisfy the necessary conditions for a
% minimum mean square error quantizer.
% nLevels - Number of quantizer output levels
% FPDF - Cell array of function handles {Farea, Fmean, Fvar}
% Farea - Pointer to a function to calculate the integral of a probability
%   density function in an interval.
%                 b
%   Farea(a,b) = Int p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fmean - Pointer to a function to calculate the mean of a probability
%   density function in an interval.
%                 b
%   Fmean(a,b) = Int x p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fvar - Pointer to a function to calculate the second moment of a
%   probability density function in an interval.
%                b
%   Fvar(a,b) = Int x^2 p(x) dx,
%                a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Sym - Logical flag. If true, the quantizer is required to be symmetric
%   about the mean. For symmetric quantizers, the probability density is
%   assumed to be symmetric about its mean.
%
% Yq - nLevels quantizer output levels (in ascending order)
% Xq - nLevels-1 quantizer decision levels (in ascending order)
% MSE - Resulting mean square error
% Entropy - Resulting entropy (bits)

Farea = FPDF{1};
Fmean = FPDF{2};

% Parameters
TolP = 1e-8;    % Symmetry check tolerance

% Consistency check
if (Sym)
  Xmean = feval(Fmean, -Inf, Inf);
  PH = feval(Farea, Xmean, Inf);
  if (abs(PH - 0.5) > TolP)
    fprintf('QuantOpt: Warning, PDF is not symmetric\n');
  end
end

% Symmetry conditions
if (~Sym)
  NlevH = nLevels;
  QSym = 0;
elseif (mod(nLevels,2) == 0)
  NlevH = nLevels/2;
  QSym = 2;
else
  NlevH = (nLevels-1)/2;
  QSym = 1;
end

% Find a minimum mean square error quantizer
Yq = QuantLloyd(NlevH, FPDF, QSym);

% Refine the solution
Yq = QuantRefine(Yq, FPDF, QSym);

% Fill in all of the levels for symmetrical distributions
if (QSym == 1)
  Yq = [(Xmean-fliplr(Yq)) Xmean Yq];
elseif (QSym == 2)
  Yq = [(Xmean-fliplr(Yq)) Yq];
end

% Generate the quantizer decision levels
if (nLevels > 1)
  Xq = 0.5 * (Yq(1:end-1) + Yq(2:end)); 
else
  Xq = [];
end

% Calculate the MSE and entropy
if (nargout > 2)
  [MSE, Entropy, SNRdB] = QuantSNR(Yq, Xq, FPDF);
end
  
return

function [Yq, Nzero] = QuantUOpt (nLevels, FPDF, Sym)
% Invoke a general search procedure to find the minimum
% mean square error uniform quantizer.
%
% The quantizer is defined by its output levels. The decision levels lie
% midway between output levels. Since the decision levels are equally spaced,
% the output levels are entirely defined by the midpoint of the quantizer and
% and quantizer step size.
%
% nLevels - Number of quantizer output levels. For symmetric quantizers
%   this is the number of levels above the mean.
% FPDF - Cell array of function pointers {Farea, Fmean, Fvar}
% Sym - Symmetry flag (optional, default false)
%
% Yq - nLevels output levels in ascending order
% Nzero - Number of intervals with zero probability

% Parameters
MaxIter = 6;

if (nargin < 3)
  Sym = false;
end

Farea = FPDF{1};
Fmean = FPDF{2};
Fvar = FPDF{3};

% Set the optimization parameters
Xmean = feval(Fmean, -Inf, Inf);
sd = feval(Fvar, -Inf, Inf) - Xmean^2;

% Convergence criteria
TolFun = sd^2 * 1e-6;
TolX = sd * 1e-6;
Options = optimset('TolFun', TolFun, 'TolX', TolX);

% Initialization
% The symmetric case is recognized in UnifQ by whether it is called
% by a one or a two element input.
V0(1) = 2 * sd / nLevels;
if (~Sym)
  V0(2) = Xmean;
end

for Iter = 1:MaxIter

% Search
  V = fminsearch(@UnifQ, V0, Options, nLevels, FPDF);

% Generate a set of uniform output values
  dY = V(1);
  if (~Sym)
    Yctr = V(2);
  else
    Yctr = Xmean;
  end
  [Yq, Xq] = GenUQuant(nLevels, dY, Yctr);

% Test for zero probability intervals
  Nzero = NzeroProb(Xq, Farea);
  if (Nzero > 0 || Iter > 1)
    fprintf('QuantUOpt: %d/%d intervals have zero probability\n', ...
            Nzero, nLevels);
  end
  if (Nzero == 0)
    break
  end

% Increase the step size for the next iteration
% If Nzero full intervals, the actual zero interval is
%   Nzero*dY < dZ < (Nzero+2)*dY
  V0(1) = nLevels/(nLevels-Nzero-1) * dY;
  
end

return

function MSE = UnifQ(V, nLevels, FPDF)
% Function to be minimized

Fmean = FPDF{2};

NV = length(V);
dY = V(1);      % Step size
if (NV == 1)
  Yctr = feval(Fmean, -Inf, Inf);
else
  Yctr = V(2);
end

% Generate a set of uniform output values
[Yq, Xq] = GenUQuant(nLevels, dY, Yctr);

MSE = QuantMSE(Yq, Xq, FPDF);

return

function [Yq, Xq] = GenUQuant (nLevels, dY, Yctr)
% Generate a set of uniform quantizer output values
%  nLevels - Number of quantizer output levels
%  dY - Quantizer step size (optional, default 1)
%  Yctr - Mean of the quantizer output levels (optional, default 0)
%
%  Xq - Quantizer decision levels (nLevels-1 values)
%  Yq - Quantizer output levels (nLevels values)
if (nargin <= 2)
  Yctr = 0;
end
if (nargin <= 1)
  dY = 1;
end

Yq = ((0:nLevels-1) - (nLevels-1)/2) * dY + Yctr;
if (nLevels > 1)
  Xq = 0.5 * (Yq(1:end-1) + Yq(2:end));
else
  Xq =[];
end

return

function Nzero = NzeroProb (Xq, Farea)

nLevels = length(Xq)+1;
Nzero = 0;
XL = -Inf;
for i = 1:nLevels
  if (i == nLevels)
    XU = Inf;
  else
    XU = Xq(i);
  end
  Area = feval(Farea, XL, XU);
  if (Area <= 0)
    Nzero = Nzero + 1;
  end
  XL = XU;
end

return

function [MSE, Entropy, SNRdB, sdV] = QuantSNR (Yq, Xq, FPDF, ScaleF)
% Calculate the SNR (dB) for a quantizer defined by a table for a
% signal with a given probability density function.
%
% Yq - Nlev quantizer output levels (in ascending order)
% Xq - Nlev-1 quantizer decision levels (in ascending order)
% FPDF - Cell array of function handles {Farea, Fmean, Fvar}
% Farea - Pointer to a function to calculate the integral of a probability
%   density function in an interval.
%                 b
%   Farea(a,b) = Int p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fmean - Pointer to a function to calculate the mean of a probability
%   density function in an interval.
%                 b
%   Fmean(a,b) = Int x p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fvar - Pointer to a function to calculate the second moment of a
%   probability density function in an interval.
%                b
%   Fvar(a,b) = Int x^2 p(x) dx,
%                a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% ScaleF - Scaling factor for the signal (optional, default 1)

% MSE - Resulting mean square error
% Entropy - Resulting entropy (bits)
% SNRdB - Resulting signal-to-noise ratio dB
% sdV - Peak to average ratio for the quantizer. The peak level (overload)
%   V of the quantizer is defined here as the the distance from the mean to
%   the largest output level plus a half step size. Then
%     sdV = sd / V,
%   where sd is the standard deviation of the signal.

if (nargin <= 3)
  ScaleF = 1;
end

% The scale factor represents a scaling of the signal
% For convenience we scale the quantizer by 1/ScaleF and later multipy
% the computed MSE by ScaleF^2 achieve the same effect as scaling the
% signal

% Quantizers are really a combination of two operations
%  - determine an interval k, such that Xq(k-1) <= x < Xq(k)
%  - output a value Yq(k)
% Consider a scaled signal g*x, then k can be determined as the value
% which satisfies
%   Xq(k-1)/g <= x < Xq(k)/g
% The quantization error for the scaled signal is
%   e = g*x - Yq(k)
%     = g * (x - Yq(k)/g)
% The bracketed term is the error for a signal x (not scaled) applied
% to a quantizer with break points Xq(.)/g and output levels Yq(.)/g.

Nlev = length(Yq);
if (ScaleF ~=1)
  if (Nlev > 1)
    Xq = Xq / ScaleF;
  end
  Yq = Yq / ScaleF;
end

% Calculate the mean square error
MSE = ScaleF^2 * QuantMSE (Yq, Xq, FPDF);

% Calculate the entropy
if (nargout >= 2)
  Entropy = QuantEntropy (Xq, FPDF{1});
end

% Calculate the SNR in dB
if (nargout >= 3)
  Fmean = FPDF{2};
  Fvar = FPDF{3};
  Xmean = feval(Fmean, -Inf, Inf);
  Svar = feval(Fvar, -Inf, Inf) - Xmean^2;
  SNRdB = 10 * log10(ScaleF^2 * Svar / MSE);
end

% Calculate the peak to average ratio for the quantizer
if (nargout >= 4)
  if (Nlev > 1)
    VL = (Yq(1)-Xq(1)) + Yq(1);
    VH = (Yq(end)-Xq(end)) + Yq(end);
    V = max(abs(VL-Xmean), abs(VH-Xmean));
  else
    V = Inf;
  end
  sdV = sqrt(Svar) / V;
end

return

function MSE = QuantMSE (Yq, Xq, FPDF)
% Calculate the mean-square quantization error
%
% This routine calculates the mean square error for a quantizer specified
% by a table of decision levels and output levels. The quantizer assigns an
% output level to a range of input values in accordance with the following
% table,
%       interval           input range          output level
%          1            -Inf < x =< Xq(1)          Yq(1)
%          2           Xq(1) < x =< Xq(2)          Yq(2)
%         ...                 ...                   ...
%        Nlev-1   Xq(Nlev-1) < x =< Xq(Nlev-1)     Yq(Nlev-1)
%        Nlev       Xq(Nlev) < x =< Inf            Yq(Nlev)
%
% Yq - Nlev quantizer output levels (in ascending order)
% Xq - Nlev-1 quantizer decision levels (in ascending order)
% FPDF - Cell array of function handles {Farea, Fmean, Fvar}
% Farea - Pointer to a function to calculate the integral of a probability
%   density function in an interval.
%                 b
%   Farea(a,b) = Int p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fmean - Pointer to a function to calculate the mean of a probability
%   density function in an interval.
%                 b
%   Fmean(a,b) = Int x p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
% Fvar - Pointer to a function to calculate the second moment of a
%   probability density function in an interval.
%                b
%   Fvar(a,b) = Int x^2 p(x) dx,
%                a
%   where (a,b) is the interval and p(x) is the probability density
%   function.
%
% The mean square error is calculated as
%         Nlev Xq(i)
%   MSE = Sum   Int  (x-Yq(i))^2 p(x) dx
%         i=1  Xq(i-1)
%
%                     Nlev          Xq(i)
%       = var[p(x)] - Sum  [2 Yq(i)  Int x p(x) dx
%                     i=1           Xq(i-1)
%
%                                     Xq(i)
%                           - Yq(i)^2  Int p(x) dx ],
%                                     Xq(i-1)
% with Xq(0)=-Inf and Xq(Nlev+1)=Inf.

Farea = FPDF{1};
Fmean = FPDF{2};
Fvar = FPDF{3};

% Sum the mean square error terms, interval by interval
Nlev = length(Yq);
if (Nlev ~= length(Xq)+1)
  error('QuantMSE - Invalid length for Yq and/or Xq');
end

dMSE = 0;
XU = -Inf;
for i = 1:Nlev
  XL = XU;
  if (i < Nlev)
    XU = Xq(i);
  else
    XU = Inf;
  end
  dMSE = dMSE + Yq(i)*(2*feval(Fmean, XL, XU) - Yq(i)*feval(Farea, XL, XU));
end

MSE = feval(Fvar, -Inf, Inf) - dMSE;

return

function Ent = QuantEntropy (Xq, Farea)
% Calculate the entropy for a quantizer
%
% This routine calculates the entropy for a quantizer specified by
% a table of decision levels. The quantizer assigns an output level
% to a range of input values in accordance with the following table,
%       interval            input range
%          1             -Inf < x =< Xq(1)
%          2            Xq(1) < x =< Xq(2)
%         ...                  ...
%        Nlev-1    Xq(Nlev-1) < x =< Xq(Nlev-1)
%        Nlev        Xq(Nlev) < x =< Inf
%
% The entropy is calculated as
%           Nlev
%   Ent = - Sum  p(i) Log2(p(i)) ,
%           i=1
%              Xq(i)
% where p(i) =  INT  p(x) dx ,
%              Xq(i-1)
% where Xq(0)=-Inf and Xq(Nlev)=Inf.
%
% Xq - Nlev-1 quantizer decision levels (in ascending order)
% Farea - Pointer to a function to calculate the integral of a
%   probability density function in an interval.
%                 b
%   Farea(a,b) = Int p(x) dx,
%                 a
%   where (a,b) is the interval and p(x) is the probability density
%   function.

% Sum the entropy, interval by interval
Nlev = length(Xq) + 1;
Enat = 0;
XU = -Inf;
for i = 1:Nlev

% Calculate the probability of the interval (XQL,XQU)
  XL = XU;
  if (i < Nlev)
    XU = Xq(i);
  else
    XU = Inf;
  end
  Pr = feval(Farea, XL, XU);

% Calculate the contribution to the entropy
  if (Pr > 0)
    Enat = Enat - Pr * log(Pr);
  elseif (Pr < 0)
    error('QuantEntropy: Negative interval probability');
  end
end

% Return the entropy in bits
Ent= Enat / log(2);

return

function Yq = QuantLloyd (Nlev, FPDF, QSym)
%  Iterate to find the output levels for a minimum mean square error
%  quantizer.
%
% This subroutine searches for a set of quantizer output levels which are
% the conditional means (centroids) of the quantizer decision regions. The
% decision boundaries are assumed to lie mid-way between output levels.
%
% The procedure is based on a Lloyd-Max iteration.
% (1) Given the start of an interval and an output level, find the end of
%     the interval, such that the output level is the centroid of the
%     probability density function in that interval.
% (2) Use the output level and upper interval edge to find the output
%     level in the next interval (the interval edge must lie midway
%     between output levels).
% (3) With the start and output level of the next interval
%     determined, repeat step (1) for this interval.
% The iteration continues by modifying the initial output level based
% on whether the centroid property of the last interval was satisfied.
%    S. P. Lloyd, "Least Squares Quantization in PCM", IEEE Trans. Inform.
%      Theory, vol. 28, no. 2, pp. 129-137, March 1982.
%    J. Max, "Quantizing for Minimum Distortion", IEEE Trans. Inform.
%      Theory, vol. 6, no. l, pp. 7-12, March 1960.
%
% Nlev - Number of quantizer output levels. For symmetric quantizers
%   this is the number of levels above the mean.
% FPDF - Cell array of function pointers {Farea, Fmean, Fvar}
% QSym - Symmetry flag (optional, default 0)
%   0 - Quantizer not constrained to be symmetric
%   1 - Quantizer is symmetric with an odd number of levels. The middle
%       level is fixed at the mean. This routine finds the Nlev output
%       levels above the mean.
%   2 - Quantizer is symmetric with an even number of levels. The middle
%       decision level is at the mean. This routine finds the Nlev output
%       levels above the mean.
%
% Yq - Nlev output levels in ascending order

% There are 3 symmetry cases to consider for symmetry
% QSym == 0: No symmetry. The first interval is defined by a lower
%   boundary XL = -Inf and an output level (centroid) Yc. The centroid
%   position is iterated.
% QSym == 1: Symmetrical quantizer with fixed output level at the mean.
%   If we place another output level Yc, the decision level XL lies
%   midway between the mean and Yc. When we iterate Yc, the decision level
%   is also be readjusted.
% QSym == 2: Symmetric quantizer with a decision boundary XL at the mean.
%   The output level Yc is iterated. 

if (nargin < 3)
  QSym = 0;
end

Fmean = FPDF{2};
Fvar = FPDF{3};

% Parameters
MaxIter = 100;
TolR = 1e-5;

% Set the optimization parameters
Xmean = feval(Fmean, -Inf, Inf);
sd = feval(Fvar, -Inf, Inf) - Xmean^2;

% Convergence criterion
Tol = TolR * sd;

if (QSym == 0)                  % No symmetry
  XL = -Inf;
  Yc = Xmean - 4 * sd;
  Xstep = sd;
elseif (QSym == 1)              % Symmetric, odd number of coefficients
  Xstep = 2 * sd / Nlev;
  Yc = Xmean + Xstep;
  XL = 0.5 * (Xmean + Yc);
elseif (QSym == 2)              % Symmetric, even number of coefficients
  Xstep = 2 * sd / Nlev;
  XL = Xmean;
  Yc = Xmean + Xstep;
end

% Initialization
Ytrial = Yc;
Ybase = Yc;
Xstep = 0.5 * Xstep;
FL = false;         % True if an upper bound has been found
FU = false;         % True if a lower bound has been found

for Iter = 1:MaxIter

% Find a set of output levels satisfying the necessary conditions
% for a minimum mean square error quantizer
% XU is the largest quantizer decision level
%  Yc = Ytrial;
  if (QSym == 1)
    XL = 0.5 * (Xmean + Ytrial);
  end

  [Yq, XU] = QuantLevel(Nlev, FPDF, XL, Ytrial);
  
  if (isnan(XU) || Yq(end) == XU)
    FU = true;
  else
    FL = true;
    Ybase = Ytrial;
  end

% Check for convergence, adjust the step size
  if (FL && FU)
    if (Xstep < Tol)
      break
    end
    Xstep = 0.5 * Xstep;
    Ytrial = Ybase + Xstep;
  elseif (FL)     % Need to extend the search upward
    Xstep = 2 * Xstep;
    Ytrial = Ybase + Xstep;
  else            % Need to back up the initial value
    Ybase = Ytrial;
    Xstep = 2 * Xstep;
    if (~isinf(XL))
      Xstep = min(Xstep, 0.5*(Ybase - XL)); % Don't step below XL
    end
    Ytrial = Ybase - Xstep;
  end

end

if (Iter >= MaxIter)
  error('QuantLloyd: Failed to converge');
end
if (~FU || ~FL)
  error('QuantLloyd: Feasible solution not found');
end
fprintf('QuantLloyd: Converged, %d iterations\n', Iter);

return

function [Yq, XU] = QuantLevel(Nlev, FPDF, XL, Yc)
% Find a set of quantizer output levels, given the lower boundary and
% centroid of the first interval.
%
% Given the lower boundary (decision level) for an initial interval and the
% output level (centroid) of that interval, this routine first finds the
% upper decision boundary for that interval. These values are telescoped to
% give the lower boundary and output level for the second interval. This
% process continues for each interval.
%
% The last decision level is returned by this routine. If this value is
% finite, the last interval does not fully encompass the tail of the
% probability density function. If the last decision level is NaN, the last
% output level is not the centroid of the last region extenting to infinity.
%
% Nlev - Number of quantizer output levels. For symmetric quantizers this
%   is the number of levels above the mean.
% FPDF - Cell array of function handles {Farea, Fmean, Fvar}
% XL - Lower boundary of the first interval
% Yc - Centroid of the first interval
%
% Yq - Nlev quantizer output levels in ascending order. Trailing values may
%   be NaN if it is not possible to have intervals with these output levels
%   as the centroids of the corresponding intervals. The quantizer decision
%   levels lie midway between output levels. 
% XU - Largest decision level (greater than or equal to Yq(Nlev)) or
%   NaN.

% If XU is finite, then the levels should be adjusted upward (increase
% Yc). If XU is NaN, the last output level is not the centroid of the last
% interval, and the levels should be adjusted downward. Another case occurs
% if an entire interval is zero probability. Then the upper decision level
% falls on the output level for that interval. Test for XU == Yq(Nlev).

Yq = zeros(1, Nlev);
for i = 1:Nlev

% Given the lower decision level and the output level for an interval,
% find the upper decision level
  if (isnan(XL))
    XU = XL;
    Yc = XL;
  else
    XU = QuantInterval(FPDF, XL, Yc);
  end

% Given the interval limits just found, telescope to find the
% output level to be used in the next iteration
  Yq(i) = Yc;
  XL = XU;
  Yc = XU + (XU - Yc);

end

return

function Xb = QuantInterval (FPDF, Xa, Yc)
% Find the upper edge of an interval which has a given centroid.
%
% This routine solves for upper limit of the integral
%    Xb
%   Int (x-Yc) p(x) dx = 0.
%    Xa
%
% The routine is designed to allow for Yc to be less than Xa, in which case
% Xb <= Yc, or for Yc to be greater than Xa, in which case Xb >= Yc.
%
% FPDF - Cell array of function handles {Farea, Fmean, Fvar}
% Xa - Lower boundary of the interval
% Yc - Centroid of the interval
%
% Xb - Returned value representing the upper boundary of the interval. This
%   value is set to NaN if there is no solution for Xb on the opposite side
%   of Yc from Xa.

% Parameters
XstepR = 1.2;    % Initial relative step size
TolF = 1e-5;     % Function amplitude relative tolerance
TolX = 1e-5;     % Position relative tolerance
MaxIter = 100;
EpsdF = 0.01;    % Choose between bisection or linear interpolation
EpsdX = 0.05;    % For linear interpolation, constrain the relative
                 % step size to EpsdX <= dX <= 1-EpsdX.

% Searching for a solution of the equation F(Xb) = 0. Consider the
% case that Yc > Xa. The function F(Xb) is zero at Xb = Xa. It becomes
% negative with increasing Xb (since p(x) is positive). It takes on its
% most negative value at Xb = Yc. It then decreases as Xb increases. It
% is the second zero crossing we seek.

% The search procedure has as its stopping criteria:
%  a) abs(F(Xb) < TolF*abs(F(YC)).
%  b) The interval of uncertainty is less than TolX*abs(Xb)
%  c) The probability from the present trial point to Inf or -Inf (as
%     appropriate) is zero. In this case Xb is set to NaN.
%  d) The maximum number of iterations is exceeded.
%  e) F(Yc)=0. This indicates that the probability density function is
%     zero in the interval (Xa,Yc).

if (Yc == Xa)
  Xb = Yc;
  return
end

Farea = FPDF{1};
Fmean = FPDF{2};
Fvar = FPDF{3};

% Evaluate the function at Yc to check if the probability
% is zero, and to determine the stopping criterion Feps
Fm = feval(Fmean, Xa, Yc) - Yc*feval(Farea, Xa, Yc);  % Should be negative
if (Fm >= 0)
  if (Fm > 0)
    error('QuantInterval: Error, function value at centroid positive');
  end
  Xb = Yc;          % Fm == 0;
  return
end

% Set up the boundaries of the search and the step size
FL = Fm;            % Value at lower boundary (initially negative)
XL = Yc;            % Lower boundary of search
if (~isinf(Xa))

  Xb = Yc + XstepR * (Yc - Xa);  % Initial trial upper boundary

else

% Find the standard deviation of the distribution
  Xmean = feval(Fmean, -Inf, Inf);
  sd = sqrt(feval(Fvar, -Inf, Inf) - Xmean^2);
  Xb = Yc + sign(Yc-Xa) * sd;   % Initial test upper boundary
end

FR = -1;            % Same sign as FL to indicate no zero found
XR = Xb;
Feps = TolF * abs(Fm); % Tolerance on integral value
Xstep = 0.5 * (Xb - Yc);

% Search loop
for Iter = 1:MaxIter

  Fp = Fm;          % Previous Fm
  Fm = feval(Fmean, Xa, Xb) - Yc*feval(Farea, Xa, Xb);

% Update the end-points of the search
if (Fm >= 0)
  FR = Fm;
  XR = Xb;
else
  FL = Fm;
  XL = Xb;
end

% ----- ------
  if (FR >= 0)
% Straddling a root

% Check for convergence in the function value
% Check for convergence of the position
    if (abs(Fm) <= Feps) || ...
       (~isinf(XL) && ~isinf(XR) && ...
        abs(XR-XL) <= TolX*max(abs(XR), abs(XL)))
      return
    end

% The root location is sought by cautious linear interpolation; however
% if successive function values are nearly equal, bisection is used.
    if (abs(Fp-Fm) > EpsdF * abs(FL-FR))
      dX = FL / (FL-FR);    % Estimated zero crossing position
      dX = max(min(dX, 1-EpsdX), EpsdX); % Avoid region close to XL or XR
    else
      dX = 0.5;             % Bisection
    end
    Xb = XL + dX*(XR-XL);

% ------ ------
  else
% Not straddling a root

% Right boundary undefined
% Check for successive identical returned values
    if (Yc > Xa)
      if (Fm == Fp && feval(Farea, Xb, Inf) <= 0)
        Xb = NaN;
        return
      end
    else
      if (Fm == Fp && feval(Farea, -Inf, Xb) <= 0)
        Xb = NaN;
        return
      end
    end

% Increase the step size
    Xb = Xb + Xstep;
    Xstep = 2 * Xstep;

  end

end

error('QuantInterval: Failed to converge');

function Yq = QuantRefine (Yq, FPDF, QSym)
% Iterate to find the output levels for a MMSE quantizer.
%
% This subroutine searches for a set of quantizer output levels which are
% the conditional means (centroids) of the quantizer decision regions. The
% decision boundaries are assumed to lie mid-way between output levels.
%
% Given a set of initial quantizer output levels, the mean square error can
% be minimized by choosing the decision levels to be midway between output
% levels. Given a set of decision levels, the MSE can be minimized by
% choosing the output levels to be the centroids of the intervals between
% decision levels. The iteration involves alternately adjusting the
% decisions levels and the output levels. At each step, the MSE will be
% non-increasing. This is the Method I iteration proposed by Lloyd.
%   S. P. Lloyd, "Least Squares Quantization in PCM", IEEE Trans. Inform.
%     Theory, vol. 28, no. 2, pp. 129-137, March 1982.

% Yq - Initial Nlev quantizer output levels in ascending order
% FPDF - Cell array of function pointers {Farea, Fmean, Fvar}
% QSym - Symmetry flag, optional (default 0)
%   0 - Quantizer not constrained to be symmetric
%   1 - Quantizer is symmetric with an odd number of levels. The middle
%       level is fixed at the mean. This routine finds the Nlev output
%       levels above the mean.
%   2 - Quantizer is symmetric with an even number of levels. The middle
%       decision level is at the mean. This routine finds the Nlev output
%       levels above the mean.
%
% Yq - Final Nlev quantizer output levels in ascending order

if (nargin < 3)
  QSym = 0;
end

Farea = FPDF{1};
Fmean = FPDF{2};
Fvar = FPDF{3};

% Parameters
MaxIter = 4000;
TolR = 1e-6;

% Find the standard deviation
Xmean = feval(Fmean, -Inf, Inf);
sd = sqrt(feval(Fvar, -Inf, Inf) - Xmean^2);

dMSE = 0;
Nzero1 = 0;
dXTol = TolR * sd;
Xmean = feval(Fmean, -Inf, Inf);
Nlev = length(Yq);

for Iter = 1:MaxIter

% Initialization for the first interval
  dMSEp = dMSE;
  Nzero = 0;
  dX = 0;
  dMSE = 0;

  if (QSym == 0)
    XL = -Inf;
  elseif (QSym == 1)
    XL = 0.5 * (Xmean + Yq(1));
% Add MSE contribution of the half interval from the mean to
% the first decision level
    Yc = Xmean;     % Output level at the mean
    Area = feval(Farea, Yc, XL);
    Avg = feval(Fmean, Yc, XL);
    dMSE = Yc * (2 * Avg - Yc * Area);
  elseif (QSym == 2)
    XL = Xmean;
  end

% -----
  for i = 1:Nlev

% Find the upper decision level
    if (i < Nlev)
      XU = 0.5 * (Yq(i+1) + Yq(i));
    else
      XU = Inf;
    end
    Area = feval(Farea, XL, XU);

% If the probability of the interval is zero, the output level is
% placed at the edge of the interval nearest the mean to move the
% output level towards a region with non-zero probability
    if (Area < 0)
      error('QuantRefine: Negative interval area');

    elseif (Area == 0)

      Nzero = Nzero + 1;
      if (abs(XL - Xmean) < abs(XU - Xmean))
        Yc = XL;
      else
        Yc = XU;
      end

    else

% Choose the output level to be the centroid of the interval
      Avg = feval(Fmean, XL, XU);
      Yc = Avg / Area;

% Accumulate the decrease in mean square error (relative to a one level
% quantizer placed at the mean) due to each quantizer interval. The MSE can
% be expressed as var(p(x))-dMSE. For symmetric quantizers, dMSE is the
% decrease in MSE for the intervals above the mean (the full dMSE is twice
% the value calculated).
      dMSE = dMSE + Yc * (2 * Avg - Yc * Area);
    end

% Find the maximum (relative / absolute) change in this iteration
    dX = max(dX, abs(Yc - Yq(i)) / max(1, abs(Yq(i))));

% Determine the next decision level and set the output level
% There are several ways to go here.
% (1) Adjust all output levels using the decision levels based on the
%     original output levels. The decision levels will get updated on the
%     next iteration cycle.
% (2) As an output level gets set, modify the decision level at the upper
%     end. This will affect the output level of the next interval in the
%     same iteration cycle.
% (Non-exhaustive) tests seem to show method (1) is best for non-symmetric
% quantizers where we move from outer levels to inner levels to outer
% levels, while method (2) is best for symmetic quantizers where we move
% from inner levels to outer levels. We opt for method (1) because we can
% calculate dMSE with essentially no additional cost and use this value as
% a safety net to halt the iterations when dMSE no longer changes.
%   XL = 0.5 * (Yq(i+1) + Yq(i)); % Method (1)
%   XL = 0.5 * (Yq(i+1) + Yc);    % Method (2)
    XL = XU;
    Yq(i) = Yc;

  end
% -----

  if (Iter == 1)
    Nzero1 = Nzero;
  end

% Check for convergence
  if (Nzero == 0 && (dX < dXTol || dMSE <= dMSEp))
    break
  end

end

if (Iter == MaxIter)
  fprintf('QuantRefine: Maximum iteration count reached\n');
end

if (dX < dXTol)
  fprintf('QuantRefine: Converged, %d iterations:', Iter);
  fprintf(' level changes less than threshold\n');
end
if (dMSE <= dMSEp)
  fprintf('QuantRefine: Converged, %d iterations:', Iter);
  fprintf(' MSE not decreasing\n');
end
if (Nzero1 > 0)
  fprintf('QuantRefine: %d zero probability interval(s) initially\n', Nzero1);
  fprintf('QuantRefine: %d zero probability interval(s) remain\n', Nzero);
end

return
