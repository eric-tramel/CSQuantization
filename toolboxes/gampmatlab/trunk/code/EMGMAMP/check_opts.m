%***********NOTE************
% The user does not need to call this function directly.  It is called
% insead by EMGMAMP.
%
% This function checks if the GAMP and EM options are set, and if not, set
% them to defaults.
%
%
% Coded by: Jeremy Vila, The Ohio State Univ.
% E-mail: vilaj@ece.osu.edu
% Last change: 4/4/12
% Change summary: 
%   v 1.0 (JV)- First release
%
% Version 2.0
function [EMopt] = check_opts(EMopt)

%Verify that EM max iter is set
if ~isfield(EMopt,'maxEMiter')
    EMopt.maxEMiter = 20;
end

%Verify that EM tolerance is set
if ~isfield(EMopt,'EMtol')
    EMopt.EMtol = 1e-5;
end;

% Default to learning active mean
if ~isfield(EMopt,'learn_mean')
    EMopt.learn_mean = true;
end;

% Default to learning active variance
if ~isfield(EMopt,'learn_var')
    EMopt.learn_var = true;
end;

if ~isfield(EMopt,'learn_weights')
    EMopt.learn_weights = true;
end;

% Default to learning lambda
if ~isfield(EMopt,'learn_lambda')
    EMopt.learn_lambda = true;
end;

% Default to learning noise variance
if ~isfield(EMopt,'learn_noisevar')
    EMopt.learn_noisevar = true;
end;

% Default to learning noise variance
if ~isfield(EMopt,'sig_dim')
    EMopt.sig_dim = 'col';
end;

%Verify that Lmax is set
if ~isfield(EMopt,'L')
    EMopt.L = 3;
end;

%Specify a defaul mode
if ~isfield(EMopt,'heavy_tailed')
    EMopt.heavy_tailed = false;
end;

return