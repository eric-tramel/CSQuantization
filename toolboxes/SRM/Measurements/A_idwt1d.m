% Written by Thong Do, JHU, 2007
% modified, Feb. 2008

function b = A_idwt1d(A,alpha,f0,f1,L,m,n)

% alpha: sparse vector which is taken from Wavelet transform of image x
% A: random projection matrix K x N
% f0, f1: lowpass and highpass filter at the synthesis part
% L: level of decomposition
% m, n: size of image
% per : take per percent of wavelet coefficients
% Return b: vector K x 1

% converting wavelet coefficients (sparse representation) into samples
u = reshape(alpha,m,n);    
im = idwt2d(u,f0,f1,L);

% converting samples into measurements
x = im(:);    
if ~isa(A, 'function_handle')
    b = A*x;
else
    b = A(x);
end
    
   