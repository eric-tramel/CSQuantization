function rx=idwt2d(x,f0,f1,maxlevel);
% rx=idwt2d(x,f0,f1,maxlevel);
% Inverse 2D Discrete Wavelet Transform

if nargin<4
	maxlevel=1
end;

[m,n]=size(x);
rx=zeros(m,n);
m=m/(2^maxlevel);n=n/(2^maxlevel);
for i=1:maxlevel
	xl=lphrec(x(1:m,1:n)',x(1:m,n+1:2*n)',f0,f1);
	xh=lphrec(x(m+1:2*m,1:n)',x(m+1:2*m,n+1:2*n)',f0,f1);
	x(1:2*m,1:2*n)=lphrec(xl',xh',f0,f1);
	m=m*2;n=n*2;
end;
rx=x;
