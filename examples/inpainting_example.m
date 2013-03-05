% Example of how to use the inpainting toolbox

clear
csq_deps('common','inpaint');

X = csq_load_data('image','lena.jpg');


block_dim = [32 32];

X1 = deblocking_filter(X,block_dim,2,0);
X2 = deblocking_filter(X,block_dim,2,1);
X3 = deblocking_filter(X,block_dim,2,2);
X4 = deblocking_filter(X,block_dim,2,3);
X5 = deblocking_filter(X,block_dim,2,4);
X6 = deblocking_filter(X,block_dim,2,5);

csq_printf('Method 1: PSNR = %0.2fdB\n',PSNR(X,X1));
csq_printf('Method 2: PSNR = %0.2fdB\n',PSNR(X,X2));
csq_printf('Method 3: PSNR = %0.2fdB\n',PSNR(X,X3));
csq_printf('Method 4: PSNR = %0.2fdB\n',PSNR(X,X4));
csq_printf('Method 5: PSNR = %0.2fdB\n',PSNR(X,X5));
csq_printf('Method 6: PSNR = %0.2fdB\n',PSNR(X,X6));
