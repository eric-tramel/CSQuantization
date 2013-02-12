%***********NOTE************
% The user does not need to call this function directly.  It is called
% insead by EMGMAMP.
%
% This function sets initializations for all unspecified GM parameters to
% defaults.
%
% Coded by: Jeremy Vila, The Ohio State Univ.
% E-mail: vilaj@ece.osu.edu
% Last change: 4/4/12
% Change summary: 
%   v 1.0 (JV)- First release
%
% Version 2.0
%%
function [lambda, omega, theta, phi, EMopt] = set_initsGM(EMopt, Y, A, M, N, T)

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

%Initialize all parameters
if ~isfield(EMopt,'active_weights');
    L = EMopt.L;
else
    L = size(EMopt.active_weights,1);
end

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

load inits.mat
%Determine signal variance of matrix or of each column
if strcmp(EMopt.sig_dim,'joint')
    temp = (norm(Y)^2/T-M*EMopt.noise_var(1,1));
    temp = repmat(temp/sum(A.multSqTr(ones(M,1)))./lambda(1,1),[N T L]);
else
    temp = (sum(abs(Y).^2,1)-M*EMopt.noise_var);
    temp = repmat(temp/sum(A.multSqTr(ones(M,1)))./lambda(1,:),[N 1 L]);
end

% Initialize Gaussian Mixture parameters
if ~EMopt.heavy_tailed
    temp2 = zeros(1,1,L);
    omega = zeros(1,1,L);
    
    %initialize active weights with pre-defined inputs or defaults
    if isfield(EMopt,'active_weights')
        if size(EMopt.active_weights,2) > 1
            omega = zeros(1,T,L);
            omega(1,:,:) = EMopt.active_weights';
        else
            omega(1,1,:) = EMopt.active_weights;
        end
    else
        omega(1,1,:) = init(L).active_weights;
    end;
    
    %initialize active variances with pre-defined inputs or defaults
    if isfield(EMopt,'active_var')
        if size(EMopt.active_var,2) > 1
            phi = zeros(1,T,L);
            phi(1,:,:) = EMopt.active_var';
        else
            phi(1,1,:) =  EMopt.active_var;
        end
    else
        temp2(1,1,:) = init(L).active_var;
        phi =  repmat(temp2,[N T 1])*12.*temp;
    end;

    %initialize active means with pre-defines inputs or defaults
    if isfield(EMopt,'active_mean')
         if size(EMopt.active_mean,2) > 1
            theta = zeros(1,T,L);
            theta(1,:,:) = EMopt.active_mean';
         else
            theta(1,1,:) = EMopt.active_mean;
         end
    else
        temp2(1,1,:) = init(L).active_mean;
        theta = repmat(temp2, [N T 1]).*sqrt(12*temp);
    end;  

%Define Heavy tailed initializations.  Override some user-defined inputs
else
    EMopt.L = 4;
    L = 4;
    omega = zeros(1,1,L);
    theta = zeros(N,T,L);
    temp2 = zeros(1,1,L);
    
    if strcmp(EMopt.sig_dim,'joint')
        temp = (norm(Y)^2/T-M*EMopt.noise_var);
        temp = repmat(temp/sum(A.multSqTr(ones(M,1)))./lambda(1,1),[N T L]);
    else
        temp = (sum(abs(Y).^2,1)-M*EMopt.noise_var);
        temp = repmat(temp/sum(A.multSqTr(ones(M,1)))./lambda(1,:),[N 1 L]);
    end
    
    temp2(1,1,:) = 1/sqrt(L)*(1:L);
    temp2 = repmat(temp2,[N,T,1]);
    phi = temp.*temp2;
    omega(1,1,:) = ones(1,L)/L;
    EMopt.learn_mean = false;
end;

%Resize all initializations to matrix form for scalar multiplications later
lambda = resize(lambda,N,T,1);
omega = resize(omega,N,T,L);
theta = resize(theta,N,T,L);
phi = resize(phi,N,T,L);

if size(EMopt.noise_var,2) == 1
    EMopt.noise_var = repmat(EMopt.noise_var,[M T]);
else
    EMopt.noise_var = repmat(EMopt.noise_var,[M 1]);
end

return