% This code demonstrates the usage of the BIHT code for a trivial one
% dimensional signal recovery. This script additionally demonstrates the
% usage of the csq_ glue functions.

% Set dependencies
clear
csq_deps('common','biht');

% Experiment variables
N = 1024;
M = 2048;
K = 20;

% Set Parameters
params.N = 1024;
params.threshold.id = 'top';
params.threshold.k = K;
params.projection.id = 'gaussian';
params.projection.subrate = M/N;
params.biht.htol = 0;
params.biht.maxIter = 3000;
params.verbose = 1;

% Generate signal
x = zeros(N,1);
p = randperm(N);
x(p(1:K)) = randn(K,1);
x = x ./ norm(x);

% Generate projection and transform
[A AT] = csq_generate_projection(params);
psi = eye(N); 
invpsi = psi;

% Unify projection and transform.
[A AT] = csq_unify_projection(A,AT,psi,invpsi);

% Generate thresholding function
T = csq_generate_threshold(params);

% Measurement
y = sign(A(x));

% BIHT Recovery
tic
xhat = biht_1d(y,A,AT,psi,invpsi,T,[],params);
biht_time = toc;

% Evaluation
mse = sum(norm(xhat-x).^2)./N;

fprintf('Recovered x in %0.2f sec. MSE = %f\n',biht_time,mse);

figure(1); cla;
subplot(2,1,1);
	hold on;
	plot(x,'-r','LineWidth',2,'DisplayName','Original');
	plot(xhat,'-b','DisplayName','Recovered');
	hold off;
	axis tight;
subplot(2,1,2);
	plot(abs(x-xhat));
	axis tight;
	ylabel('|x - xhat|');




