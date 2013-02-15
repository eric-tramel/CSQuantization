function dx=dwt1d(x,h0,h1,maxlevel);
% dx=dwt1d(x,h0,h1,level);
% Forward 1D Discrete Symmetric Biorthogonal Wavelet Transform

if nargin<4
	maxlevel=1;
end;

n=length(x);
for i=1:maxlevel
	[xl,xh]=lphdec(x(1:n),h0,h1);
	x(1:n)=[xl xh];
	n=n/2;
end;
dx=x;
