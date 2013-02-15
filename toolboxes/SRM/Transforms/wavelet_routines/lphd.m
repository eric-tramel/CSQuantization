% LPHD         Two-dimensional wavelet transform using linear-phase biorthogonal
%              two-channel filter banks (trivial convolution method)
%
%              Usage:
%                 y=lphd(x,h0,h1,lnumber{,extmode});
%
%              x is wavelet transformed to y with lnumber iterations using
%              the ANALYSIS filter pair {h0,h1} (normalized to 1/sqrt(2)). 
%              The dimensions of x must be multiples of 2^lnumber. 
%              extmode specifies the extension mode:
%                 0 - zero padding (default)
%                 1 - symmetric extension
%                 2 - circular convolution
%
%              See also LPHR.

error ('MEX-file not accessible or not in first place of search hierarchy!');

