function mydisplay(x,xmax);
% mydisplay(x,xmax);
% Function to display an image x with pixel value range from 0 to xmax

[m,n]=size(x);
if (nargin==1)
   xmax=max(max(x));
end;
figure;
whitebg('w');
image(x);
daspect([1 1 1]);
axis('off');
if size(xmax,2)==3
   colormap(xmax);
else
   colormap(gray(xmax));
end;
