function [xl,xh]=lphdec(x,h0,h1,mode);
% [xl,xh] = lphdec(x,h0,h1,mode)
% Decomposition of given scaling coefficients into their scaling and wavelet 
% parts using a LP PR analysis bank
%	Input:
%		x   : input scaling coeffs. (vector [1,x_len] or 
%		      matrix [m,x_len], i.e., list of m vectors [1,x_len])
%		h0  : analysis filter coeffs. for H0(z) (vector [1,h0_len])
%		h1  : analysis filter coeffs. for H1(z) (vector [1,h1_len])
%		mode: Extension Mode:
%			0 - zero extension 
%			1 - symmetric extension
%			2 - circular convolution
%	Output:
%		xl  : scaling coeffs. at next coarser level (low pass)
%		      (vector [1,x_len/2] or matrix [m,x_len/2])
%		xh  : wavelet coeffs. at next coarser level (high pass)
%		      (vector [1,x_len/2] or matrix [m,x_len/2])

if nargin<4
	mode=1;	
end;

[m,x_len]=size(x);
h0_len=length(h0);
h1_len=length(h1);
if rem(h0_len,2)~=rem(h1_len,2)
	error('error: filter lengths must be EE or OO!');
end;

ext=fix(max(h0_len,h1_len)/2);	% extension size
tb=rem(h0_len,2);	% change extension type for EE- or OO-FB
if mode==1
	x=[fliplr(x(:,1+tb:ext+tb)),x,fliplr(x(:,x_len-ext+1-tb:x_len-tb))];
elseif mode==2
	x=[x(:,x_len-ext+1:x_len),x,x(:,1:ext)];
else
	x=[zeros(m,ext),x,zeros(m,ext)];
end;

len=fix((x_len+1)/2);
xh=zeros(m,len);
xl=zeros(m,len);
s=(h0_len-h1_len)/2;
if s>=0
	k1=2-tb;k2=2+s;
else
	k2=2;k1=2-tb-s;
end;
h0_ord=h0_len-1;h1_ord=h1_len-1;
for i=1:len
	xl(:,i)=x(:,k1:k1+h0_ord)*h0';
	xh(:,i)=x(:,k2:k2+h1_ord)*h1';
	k1=k1+2;
	k2=k2+2;
end;

