function [A AT] = csq_unify_projection(inPhi,inPhiT,inPsi,inPsiT)
% A = unify_projection(Phi,Psi)
% This is a convenience function which joins together the projection matrix 
% as well as the sparse transform into a single pair of function handles.
%
% Inputs:
%   Phi         --  A MxN matrix or function handle mapping N->M.
%   PhiT        --  A NxM matrix or function handle mapping M->N.
%   Psi         --  A NxD matrix (D>N) or function handle mapping N->D.
%   PsiT        --  A DxN matrix (D>N) or function handle mapping D->N.
% 
% Outputs:
%   A           -- A function handle mapping D->M.
%   AT          -- A function handle mapping M->D.

%% Input Handling
if isa(inPhi,'function_handle')
    Phi = inPhi;
else
    Phi = @(x) inPhi*x;
end

if isa(inPhiT,'function_handle')
    PhiT = inPhiT;
else
    PhiT = @(x) inPhiT*x;
end

if isa(inPsi,'function_handle')
    Psi = inPsi;
else
    Psi = @(x) inPsi*x;
end

if isa(inPsiT,'function_handle')
    PsiT = inPsiT;
else
    PsiT = @(x) inPsiT*x;
end

%% Set Handles
A =  @(x) Phi(PsiT(x));
AT = @(x) Psi(PhiT(x));