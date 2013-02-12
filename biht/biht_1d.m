function x = biht_1d(y,Ain,params)
% biht_1d(y,Phi,params)
% Code to implement the basic 1D functionality of the the
% 1-Bit CS recovery procedure: Binary Iterated Hard Thresholding.
%
% Inputs:
%	y		--		A Mx1 vector of 1-bit quantized CS measurements
%       A		--		A MxN projection matrix used to calculate y = Ax
%       params	-- 		Structure lgorithm parameters:
%					params.htol
%					params.k
%					params.ATrans
%					params.maxIter
%
% This code is based on the work presented in 
% L. Jauqes, J. Laska, P. T. Boufouros, and R. G. Baraniuk,
% "Robust 1-Bit Compressive Sensing via Binary Stable Embeddings"


%% Variable Initialization
% Flags
FLAG_bitsspecified = 0;
FLAG_Nspecified = 0;
FLAG_ATransspecified = 0;
FLAG_Kspecified = 0;
FLAG_htolspecified = 0;

%% Check the input parameters
if isfield(params,'htol')
	FLAG_htolspecified = 1;
end

if isfield(params,'k')
	FLAG_Kspecified = 1;
end

if isfield(params,'bits')
	FLAG_bitsspecified = 1;
end

if isfield(params,'N')
	FLAG_Nspecified = 1;
end

if isfield(params,'ATrans')
	FLAG_ATransspecified = 1;
end


%% Input handling
if isa(Ain,'function_handle')
	A =@(x) sign(Ain(x));
else
	A = @(x) sign(A*x);
end

if (FLAG_ATransspecified)
	if isa(params.ATrans,'function_handle')
		AT = @(x) params.ATrans(x);
	else
		AT = @(x) params.ATrans*x;
	end
else
	if isa(Ain,'function_handle')
		error('biht_1d:NoTranspose','No projection transpose function specified.');
	else
		AT = @(x) Ain'*x;
	end
end


