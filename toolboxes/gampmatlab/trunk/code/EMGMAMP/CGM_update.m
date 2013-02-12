%***********NOTE************
% The user does not need to call this function directly.  It is called
% insead by EMGMAMP.
%
% This function returns complex EM updates of the parameters lambda, 
% omega, theta, and phi given the GAMP outputs Rhat and Rvar
%
%
% Coded by: Jeremy Vila, The Ohio State Univ.
% E-mail: vilaj@ece.osu.edu
% Last change: 4/4/12
% Change summary: 
%   v 1.0 (JV)- First release
%   v 2.0 (JV)- Accounts for MMV model.  
%
% Version 2.0
%%
function [lambda, omega, theta, phi] = CGM_update(Rhat, Rvar, lambda, omega, theta, phi, EMopt)

%Calcualte Problem dimensions
L = size(omega,3);
[N, T] = size(Rhat);

D_l = zeros(N,T,L); a_nl = zeros(N,T,L);
gamma = zeros(N,T,L); nu = zeros(N,T,L);

%Evaluate posterior likelihoods
for i = 1:L
    dummy = Rvar+phi(:,:,i);
    D_l(:,:,i) = lambda.*omega(:,:,i)./dummy...
        .*exp(-abs(theta(:,:,i)-Rhat).^2./dummy);
    gamma(:,:,i) = (Rhat.*phi(:,:,i)+Rvar.*theta(:,:,i))./dummy;
    nu (:,:,i) = (Rvar.*phi(:,:,i))./dummy;
    a_nl(:,:,i) = Rvar./(dummy).*omega(:,:,i)...
    .*exp(-abs(Rhat-theta(:,:,i)).^2./dummy+abs(Rhat).^2./Rvar);
end;

%Find posterior that the component x(n,t) is active
if EMopt.learn_lambda
    a_n = real(lambda./(1-lambda).*sum(a_nl,3));
    a_n = 1./(1+a_n.^(-1));
    a_n(isnan(a_n)) = 0.001;
    if strcmp(EMopt.sig_dim,'joint')
        lambda = repmat(sum(sum(a_n))/N/T,[N,T]);
    else
        lambda = repmat(sum(a_n)/N,[N 1]);
    end
end

%Find the Likelihood that component n,t belongs to class l and is active
E_l = D_l./repmat((sum(D_l,3)+(1-lambda)./Rvar.*exp(-abs(Rhat).^2./Rvar)),[1 1 L]);
%Ensure real valued probability
E_l(isnan(E_l)) = 0.999;

%Update parameters based on EM equations
if strcmp(EMopt.sig_dim,'joint')
    if ~(T == 0 || T == 1 || N ==0 || N == 1)
        N_l = sum(sum(E_l));
        if EMopt.learn_mean
            theta = resize(sum(sum(E_l.*gamma))./N_l,N,T,L);
        end
        if EMopt.learn_var
            phi = resize(sum(sum(E_l.*(nu+abs(gamma-theta).^2)))./N_l,N,T,L);
        end
        if EMopt.learn_weights
            omega = N_l/N/T;
            omega = omega./repmat(sum(omega, 3), [1, 1, L]);
        end
    end
else
    if ~(N == 0 || N == 1)
        N_l = sum(E_l);
        if EMopt.learn_mean
            theta = resize(sum(E_l.*gamma)./N_l,N,T,L);
        end
        if EMopt.learn_var
            phi = resize(sum(E_l.*(nu+abs(gamma-theta).^2))./N_l,N,T,L);
        end
        if EMopt.learn_weights
            omega = N_l/N;
            omega = omega./repmat(sum(omega, 3), [1, 1, L]);
        end
    end
end

omega = resize(omega,N,T,L);

return;