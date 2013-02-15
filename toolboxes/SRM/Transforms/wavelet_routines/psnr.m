function y = psnr(im1,im2);
[m,n] = size(im1);
x1 = double(im1(:));
x2 = double(im2(:));
mse = norm(x1-x2);
mse = (mse*mse)/(m*n);
if mse >0
    y = 10*log10(255^2/mse);
else
    disp('infinite psnr');
end