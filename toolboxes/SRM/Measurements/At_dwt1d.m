% Written by Thong Do, JHU, 2007
% modified, Feb. 2008

function alpha = At_dwt1d(At,b,h0,h1,L,m,n)
% alpha: sparse vector which is taken from Wavelet transform of image x
% A: random projection matrix K x N
% h0, h1: lowpass and highpass filter at the analysis part
% L: level of decomposition
% m, n: size of image
% Return b: vector K x 1

% converting measurements into samples
if ~isa(At, 'function_handle')
    x = At*b;
else
    x = At(b);
end
% converting samples into wavelet coefficients (sparse representation)
im = reshape(x,m,n);
u = dwt2d(im,h0,h1,L);
alpha = u(:);


