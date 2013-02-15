function [xl,xh]=lpdec(x,h0,h1);
% [xl,xh] = lpdec(x,h0,h1)
% Decomposition of given scaling coefficients into their scaling and wavelet 
% parts using a LP PR analysis bank
%	Input:
%		x   : input scaling coeffs. (vector [1,x_len] or 
%		      matrix [m,x_len], i.e., list of m vectors [1,x_len])
%		h0  : filter coeffs. for H0(z) (vector [1,h0_len])
%		h1  : filter coeffs. for H1(z) (vector [1,h1_len])
% 		The filters are assumed to be aligned by padding zeros.
%	Output:
%		xl  : scaling coeffs. at next coarser level (low pass)
%		      (vector [1,x_len/2] or matrix [m,x_len/2])
%		xh  : wavelet coeffs. at next coarser level (high pass)
%		      (vector [1,x_len/2] or matrix [m,x_len/2])

[m,x_len]=size(x);
h0_len=length(h0);
h1_len=length(h1);
ext=floor(max(h0_len,h1_len)/2);	% extension size
fbt=rem(h0_len,2);	% what type of FB ? EE or OO ?

% symmetric extension
x=[fliplr(x(:,1+fbt:ext+fbt)),x,fliplr(x(:,x_len-ext+1-fbt:x_len-fbt))];

len=fix((x_len+1)/2); % subband size
xh=zeros(m,len);
xl=zeros(m,len);
s=h0_len-h1_len;
if fbt==1
  if s>=0 % h0 is longer
	 k0=1;k1=1+s;
  else % h1 is longer
	 k0=1-s;k1=1;
 end;
else
  if s>=0 % h0 is longer
	 k0=2;k1=2+s;
  else % h1 is longer
	 k0=2-s;k1=2;
 end;
end;
h0_ord=h0_len-1;h1_ord=h1_len-1;
h0=fliplr(h0); h1=fliplr(h1);
for i=1:len
	xl(:,i)=x(:,k0:k0+h0_ord)*h0';
	xh(:,i)=x(:,k1:k1+h1_ord)*h1';
	k0=k0+2;
	k1=k1+2;
end;

