%***********NOTE************
% The user does not need to call this function directly.  It is called
% insead by EMBGAMP.
%
% This function returns real EM updates of the parameters lambda, 
% theta, and phi given the GAMP outputs Rhat and Rvar
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
function [lambda, theta, phi] = BG_update(Rhat, Rvar, lambda, theta, phi, EMopt)

%Calcualte Problem dimensions
[N, T] = size(Rhat);

%Calculate posterior means and variances
temp = phi + Rvar;
gamma_n = (Rhat.*phi+Rvar.*theta)./temp;
nu_n = Rvar.*phi./temp;

%Find posterior that the component x(n,t) is active
A_n = (1-lambda)./(lambda).*sqrt(phi./nu_n).*exp(((Rhat-theta).^2./temp-Rhat.^2./Rvar)/2);
A_n = 1./(1+A_n);
A_n(isnan(A_n)) = 0.999;
temp2 = sum(A_n);

%Update BG parameters
if EMopt.learn_mean
    if strcmp(EMopt.sig_dim,'joint')
        theta = sum(sum(A_n.*gamma_n))/sum(temp2);
        theta = repmat(theta,[N T]);
    else
        theta = sum(A_n.*gamma_n)./temp2;
        theta = repmat(theta,[N 1]);
    end
end

if EMopt.learn_var
    if strcmp(EMopt.sig_dim,'joint')
        phi = sum(sum(A_n.*(nu_n+(gamma_n-theta).^2)))/sum(temp2);
        phi = repmat(phi,[N T]);
    else
        phi = sum(A_n.*(nu_n+(gamma_n-theta).^2))./temp2;
        phi = repmat(phi,[N 1]);
    end
end

if EMopt.learn_lambda
    if strcmp(EMopt.sig_dim,'joint')
        lambda = sum(temp2)/N/T;
        lambda = repmat(lambda, [N T]);
    else
        lambda = temp2/N;
        lambda = repmat(lambda, [N 1]);
    end
end

return;