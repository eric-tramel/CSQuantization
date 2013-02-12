%***********NOTE************
% The user does not need to call this function directly.  It is called
% insead by EMBGAMP.
%
% This function sets initializations for all unspecified BG parameters to
% defaults.
%
% Coded by: Jeremy Vila, The Ohio State Univ.
% E-mail: vilaj@ece.osu.edu
% Last change: 4/4/12
% Change summary: 
%   v 2.0 (JV)- First release
%
% Version 2.0

function [lambda, theta, phi, EMopt] = set_initsBG(EMopt, Y, A, M, N, T)

%Verify that noise variance is initialized
if ~isfield(EMopt,'noise_var')
    if strcmp(EMopt.sig_dim,'joint')
        EMopt.noise_var = norm(Y)^2/T/(M*101);
    else
        EMopt.noise_var = sum(abs(Y).^2,1)/(M*101);
    end
end;

%Determine undersampling ratio
del = M/N;

%Define cdf/pdf for Gaussian
normal_cdf = @(x) 1/2*(1 + erf(x/sqrt(2)));

%Define density of Gaussian
normal_pdf = @(x) 1/sqrt(2*pi)*exp(-x.^2/2);

temp = linspace(0,10,1024);
rho_SE = (1 - (2/del)*((1+temp.^2).*normal_cdf(-temp)-temp.*normal_pdf(temp)))...
    ./(1 + temp.^2 - 2*((1+temp.^2).*normal_cdf(-temp)-temp.*normal_pdf(temp)));
rho_SE = max(rho_SE);

%Initialize lambda
if isfield(EMopt,'lambda')
    lambda = EMopt.lambda;
else
    lambda = min(del*rho_SE,1);
end;

%Initialize active mean parameter
if isfield(EMopt,'active_mean')
    theta = EMopt.active_mean;
else
    theta = 0;
end

%Initialize active variance parameter
if isfield(EMopt,'active_var')
    phi = EMopt.active_var;
else
    if strcmp(EMopt.sig_dim,'joint')
        phi = (norm(Y)^2/T-M*EMopt.noise_var(1,1))/sum(A.multSqTr(ones(M,1)))./lambda;
    else
        phi = (sum(abs(Y).^2,1)-M*EMopt.noise_var)/sum(A.multSqTr(ones(M,1)))./lambda;
    end
end

%Resize all parameters to N by T matrix
[~,temp] = size(lambda);
if temp == T
    lambda = repmat(lambda, [N 1]);
else
    lambda= repmat(lambda, [N, T]);
end

[~,temp] = size(theta);
if temp == T
    theta = repmat(theta, [N 1]);
else
    theta = repmat(theta, [N, T]);
end

[~,temp] = size(phi);
if temp == T
    phi = repmat(phi, [N 1]);
else
    phi = repmat(phi, [N, T]);
end

if length(EMopt.noise_var) == T
    EMopt.noise_var = repmat(EMopt.noise_var,[M,1]);
end

return