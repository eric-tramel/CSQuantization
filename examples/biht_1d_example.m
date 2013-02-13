% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

% Set dependencies
clear
csq_deps('common','biht');

% Experiment variables
N = 1024;
M = 256;
K = 20;

% Generate signal
x = zeros(N,1);
p = randperm(N);
x(p(1:K)) = rand(K,1);
x = x ./ norm(x);

% Generate projection and transform
Phi = orth(randn(N,M))';
Psi = eye(N);       % Sparse in canonical basis
[A AT] = csq_unify_projection(Phi,Phi',Psi,Psi');

% Measurement
y = sign(A(x));

% Recovery parameters
params.k = K;
params.htol = 0;
params.maxIter = 3000;
params.ATrans = AT;
params.N = N;
params.verbose = 1;

% BIHT Recovery
tic
xhat = biht_1d(y,A,params);
biht_time = toc;

% Evaluation
mse = sum(norm(xhat-x).^2)./N;

fprintf('Recovered x in %0.2f sec. MSE = %f\n',biht_time,mse);

