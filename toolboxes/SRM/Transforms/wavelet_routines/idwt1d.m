function rx=idwt1d(x,f0,f1,maxlevel);
% rx=idwt1d(x,f0,f1,maxlevel);
% Inverse 1D Discrete Symmetric Biorthogonal Wavelet Transform

if nargin<4
	maxlevel=1
end;

n=length(x);
rx=zeros(1,n);
n=n/(2^maxlevel);
for i=1:maxlevel
	x(1:2*n)=lphrec(x(1:n),x(n+1:2*n),f0,f1);
        n=n*2;
end;
rx=x;
