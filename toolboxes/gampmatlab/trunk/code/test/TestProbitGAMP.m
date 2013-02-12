% TestProbitGAMP
%
% This script will run a single test of the performance of GAMP on a sparse
% binary classification problem using the ProbitEstimOut class.  In this
% test, it is assumed that a length-M binary vector, y (y(m) \in {0,1}), is
% associated with a group of N features, organized into an M-by-N feature
% matrix, A.  The entries of A are drawn i.i.d. Normal(0, 1/M).  The
% relationship between the feature matrix, A, and the class labels, y, is
% determined by a synthetically generated length-N "true" weight vector, x.
% Specifically, this relationship is governed by the equation
%               y = (1/2) * ( sgn(A*x + w) + 1 ),
% i.e., y(m) = 1 when A(m,:)*x + w(m) is greater than 0, otherwise y(m) =
% 0.  The weight vector, x, is a sparse vector whose elements are drawn as
% i.i.d. Bernoulli-Gaussian random variables, with standard-normal non-zero
% entries, i.e.,
%       p(x(n)) = (1 - K/N)*delta(x(n)) + (K/N)*Normal(x(n); 0, 1),
% where K < M is the number of important (discriminative) features.
% Entries of additive class label noise, w(m), are i.i.d. 
% Normal(0, sigma2).
%

%
% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu
% Last change: 05/30/12
% Change summary: 
%       - Created (05/29/12; JAZ)
% Version 0.1
%

Stream = RandStream.getDefaultStream;
SavedState = Stream.State;      % Save the current RNG state?
Stream.State = SavedState;      % Reset the RNG to previous initial state?

%% Define test parameters

% Problem setup
N = 512;            % # of features
NoverM = 4;         % Ratio of # of features-to-# of training samples
MoverK = 5;         % Ratio of # of training samples-to-# of "active" features
sigma2 = 1e-8;      % Variance of corrupting label noise
varyK = false;      % Draw K from binomial dist. (true) or keep fixed (false)
Mtest = 256;        % # of test set samples for computing error rate
    
% Probit regression function parameters
% The output channel is specified as p(y(m) = 1 | z) = Phi((z-m)/sqrt(v)),
% where Phi(*) is the CDF of a normal distribution with mean m and variance
% v.
probit_mean = 0;    % Mean m of probit regression function
probit_var = .1;     % Variance v of probit regression function


%% Generate a synthetic trial

Mtrain = round(N / NoverM);      % # of training samples
K = round(Mtrain / MoverK);      % # of "active" (discriminative) features

% First construct the "true" weight vector that yields the class labels
switch varyK
    case true
        x = zeros(N,1);
        K_true = binornd(N, K/N);   % Draw K from binomial distribution
        active_features = randperm(N, K_true);  % Indices of discriminative features
        x(active_features) = randn(K_true,1);
    case false
        x = zeros(N,1);
        K_true = K;     % Fix K_true to equal average K
        active_features = randperm(N, K_true);  % Indices of discriminative features
        x(active_features) = randn(K_true,1);
    otherwise
        error('Unrecognized option: varyK')
end

% Now construct the label noise vector, the random feature matrix, and true
% class labels for the training set
A_train = (1/sqrt(Mtrain))*randn(Mtrain,N);         % Feature matrix
w_train = sqrt(sigma2)*randn(Mtrain,1);             % Label noise
y_train = (1/2)*(sign(A_train*x + w_train) + 1);	% True binary {0,1} class labels

% Finally, construct the label noise vector, the random feature matrix, and
% true class labels for the training set
A_test = (1/sqrt(Mtrain))*randn(Mtest,N);       % Feature matrix (use training variance)
w_test = sqrt(sigma2)*randn(Mtest,1);           % Label noise
y_test = (1/2)*(sign(A_test*x + w_test) + 1);	% True binary {0,1} class labels


%% Recover with GAMP

% Build the EstimIn object for a Bernoulli-Gaussian signal
EstimIn = AwgnEstimIn(zeros(N,1), ones(N,1));
EstimIn = SparseScaEstim(EstimIn, K/N);

% Build the EstimOut object for a probit regression channel
EstimOut = ProbitEstimOut(y_train, probit_mean, probit_var);
% EstimOut = LogitEstimOut(y_train);

% Specify GAMP options
GAMPopts = GampOpt();
GAMPopts.adaptStep = false;

% Execute GAMP
[x_hat, x_hat_var] = gampEst(EstimIn, EstimOut, A_train, GAMPopts);


%% Compute error metrics

y_hat_train = (1/2) * (sign(A_train*x_hat) + 1);
y_hat_test = (1/2) * (sign(A_test*x_hat) + 1);

train_error = sum(y_hat_train ~= y_train) / Mtrain;
test_error = sum(y_hat_test ~= y_test) / Mtest;

fprintf('%% of class 1 in training set: %1.1f\n', 100*sum(y_train == 1)/Mtrain)
fprintf('%% training samples misclassified: %1.1f\n', 100*train_error)
fprintf('%% test samples misclassified: %1.1f\n', 100*test_error)


%% Plot the truth and the recovery

figure(1); clf
stem(x, 'b'); hold on; stem(x_hat, 'r'); hold off
title(sprintf(['%% training samples misclassified: %1.1f | ' ...
    '%% testing samples misclassified: %1.1f'], 100*train_error, ...
    100*test_error))
xlabel('Feature [n]'); ylabel('Feature weight')
legend('"True" weight vector', 'GAMP weight vector', 'Location', 'Best')