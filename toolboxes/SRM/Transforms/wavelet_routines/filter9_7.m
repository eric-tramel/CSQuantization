function [h0,h1,f0,f1] = filter9_7()

p0=[-5 0 49 0 -245 0 1225 2048 1225 0 -245 0 49 0 -5]/2048; %half-band filter
r=roots(p0); b0=[1 8 28 56 70 56 28 8 1];
q0=deconv(p0,b0);
r=[roots(q0); -ones(8,1)];
h0=poly([r(2:5)' -1 -1 -1 -1]);
f0=poly([r(1) r(6) -1 -1 -1 -1]);
h0=real(sqrt(2)*h0/sum(h0));
f0=real(sqrt(2)*f0/sum(f0));
f1=h0.*(-1).^[0:length(h0)-1];
h1=f0.*(-1).^[1:length(f0)];