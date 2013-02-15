function x=lprec(xl,xh,f0,f1);
% x = lprec(xl,xh,f0,f1);
% Reconstruction of given scaling and wavelet components into the higher level 
% scaling coeffs. using a LP PR synthesis bank
%	Input:
%		xl  : input scaling coeffs.
%		      (vector [1,x_len/2] or list of vectors [m,x_len/2])
%		xh  : input wavelet coeffs.
%		      (vector [1,x_len/2] or list of vectors [m,x_len/2])
%		f0  : filter coeffs. for F0(z) (vector [1,h0_len])
%		f1  : filter coeffs. for F1(z) (vector [1,h1_len])
% 		The filters are assumed to be aligned by padding zeros.
%	Output:
%		x   : scaling coeffs. at next higher level
%		      (vector [1,x_len] or list of vectors [m,x_len])

[m,xl_len]=size(xl);
f0_len=length(f0);
f1_len=length(f1);

x_len=xl_len*2;
fbt=rem(f0_len,2);
xl2=zeros(m,x_len);
xh2=zeros(m,x_len);
% upsampling
xl2(:,2:2:x_len)=xl;
xh2(:,2:2:x_len)=xh;

ext=fix(max(f0_len,f1_len)/2); % symmetric extension
if fbt==1
%xl2=[fliplr(xl2(:,2:ext+1)),xl2,zeros(m,1),fliplr(xl2(:,x_len-ext+1:x_len))];
%xh2=[fliplr(xh2(:,2:ext+1)),xh2,fliplr(xh2(:,x_len-ext:x_len-1))];
  xl2=[fliplr(xl2(:,2:ext+1)),xl2,fliplr(xl2(:,x_len-ext:x_len-1))];
  xh2=[fliplr(xh2(:,2:ext+1)),xh2,fliplr(xh2(:,x_len-ext:x_len-1))];
else
  xl2=[fliplr(xl2(:,2:ext)),xl2,zeros(m,1),fliplr(xl2(:,x_len-ext+1:x_len))];
  xh2=[-fliplr(xh2(:,2:ext)),xh2,zeros(m,1),-fliplr(xh2(:,x_len-ext+1:x_len))];
end;

xl2
xh2

x=zeros(m,x_len);
f0_ord=f0_len-1; f1_ord=f1_len-1;
s=f0_ord-f1_ord
f0=fliplr(f0); f1=fliplr(f1);
if s>=0 % if f0 is longer
	for i=1:x_len
	  x(:,i)=xl2(:,i:i+f0_ord)*f0'+xh2(:,i+s:i+s+f1_ord)*f1';
	end;
else % if f0 is shorter
	for i=1:x_len
	  x(:,i)=xl2(:,i-s:i-s+f0_ord)*f0'+xh2(:,i:i+f1_ord)*f1';
	end;
end;
