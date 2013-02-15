function dx=dwt2d(x,h0,h1,maxlevel);
% dx=dwt2d(x,h0,h1,level);
% Forward 2D Discrete Wavelet Transform

if nargin<4
	maxlevel=1;
end;

[m,n]=size(x);
for i=1:maxlevel
	[xl,xh]=lphdec(x(1:m,1:n),h0,h1);
	[xll,xlh]=lphdec(xl',h0,h1);
	[xhl,xhh]=lphdec(xh',h0,h1);
	x(1:m,1:n)=[xll xhl;xlh xhh]';
	m=m/2;n=n/2;
end;
dx=x;
