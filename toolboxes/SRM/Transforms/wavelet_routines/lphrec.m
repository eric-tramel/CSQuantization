function x=lphrec(xl,xh,f0,f1,mode);
% x = lphrec(xl,xh,f0,f1,mode)
% Reconstruction of given scaling and wavelet components into the higher level 
% scaling coeffs. using a LP PR synthesis bank
%	Input:
%		xl  : input scaling coeffs.
%		      (vector [1,x_len/2] or list of vectors [m,x_len/2])
%		xh  : input wavelet coeffs.
%		      (vector [1,x_len/2] or list of vectors [m,x_len/2])
%		f0  : synthesis filter coeffs. for F0(z) (vector [1,h0_len])
%		f1  : synthesis filter coeffs. for F1(z) (vector [1,h1_len])
%		mode: Extension Mode:
%			0 - zero extension
%			1 - symmetric extension
%			2 - circular convolution
%	Output:
%		x   : scaling coeffs. at next higher level
%		      (vector [1,x_len] or list of vectors [m,x_len])

if nargin<5
	mode=1;
end;

[m,xl_len]=size(xl);
h0_len=length(f0);
h1_len=length(f1);
if rem(h0_len,2)~=rem(h1_len,2)
	error('error: filter lengths must be EE or OO!');
end;
if [m,xl_len]~=size(xh)
	error('error: both channels must have same size!');
end;

x_len=xl_len*2;
tb=rem(h0_len,2);
ta=1-tb;
xl2=zeros(m,x_len);
xh2=zeros(m,x_len);
xl2(:,2-tb:2:x_len)=xl;
xh2(:,2:2:x_len)=xh;
ext=fix(max(h0_len,h1_len)/2);

if mode==1
	if tb==1
	  xl2=[fliplr(xl2(:,2:ext+1)),xl2,fliplr(xl2(:,x_len-ext:x_len-1))];
	  xh2=[fliplr(xh2(:,2:ext+1)),xh2,fliplr(xh2(:,x_len-ext:x_len-1))];
	else
	  xl2=[fliplr(xl2(:,2:ext)),xl2,zeros(m,1),fliplr(xl2(:,x_len-ext+1:x_len))];
	  xh2=[-fliplr(xh2(:,2:ext)),xh2,zeros(m,1),-fliplr(xh2(:,x_len-ext+1:x_len))];
	end;
elseif mode==2
	xl2=[xl2(:,x_len-ext+1+ta:x_len),xl2,xl2(:,1:ext+ta)];
	xh2=[xh2(:,x_len-ext+1+ta:x_len),xh2,xh2(:,1:ext+ta)];

else
	xl2=[zeros(m,ext-ta),xl2,zeros(m,ext+ta)];
	xh2=[zeros(m,ext-ta),xh2,zeros(m,ext+ta)];
end;

x=zeros(m,x_len);
f0_ord=h0_len-1;f1_ord=h1_len-1;
s=(f0_ord-f1_ord)/2;
if s>=0
	for i=1:x_len
	  x(:,i)=xl2(:,i:i+f0_ord)*f0'+xh2(:,i+s:i+s+f1_ord)*f1';
	end;
else
	for i=1:x_len
	  x(:,i)=xl2(:,i-s:i-s+f0_ord)*f0'+xh2(:,i:i+f1_ord)*f1';
	end;
end;
