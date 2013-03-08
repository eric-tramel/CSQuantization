clear;
csq_deps('srm','common');
subrate = 0.5;
blksize = 32;
trans_mode = 'BWHT';
block_dim = [32 32];
timing_trials = 20;


X = csq_load_data('image','lena.jpg');
x = X(:);

imsize = size(X);
N = imsize(1)*imsize(2);


% No blocks
[A AT] = projection_srmblk(subrate,N,trans_mode,blksize);

% Block-based
[ABB ATBB] = projection_srmblk(subrate,N,trans_mode,blksize,1,imsize,block_dim);


% Measurements
y  = A(x);
xt = AT(y);

yBB  = ABB(x);
xtBB = ATBB(yBB);

d_psnr   = PSNR(x,xt);
d_psnrBB = PSNR(x,xtBB);

% Display
figure(1);
subplot(1,2,1);
    imagesc(reshape(xt,imsize));
    title('Global Measurements');
    xlabel(sprintf('Back Projection\nSubrate = %0.2f\nPSNR = %0.2fdB',subrate,d_psnr));
    axis image;
subplot(1,2,2);
    imagesc(reshape(xtBB,imsize));
    title('Block Measurements');
    xlabel(sprintf('Back Projection\nSubrate = %0.2f\nPSNR = %0.2fdB',subrate,d_psnrBB));
    axis image;
colormap(gray);

% Timings
tic
for i=1:timing_trials
    y = A(x);
end
Atime = toc;

tic
for i=1:timing_trials
    xt = AT(y);
end
ATtime = toc;

tic
for i=1:timing_trials
    yBB = ABB(x);
end
ABBtime = toc;

tic
for i=1:timing_trials
    xtBB = ATBB(yBB);
end
ATBBtime = toc;

fprintf('\n-------------------\n');
fprintf('A    = %f sec\n',Atime./timing_trials);
fprintf('AT   = %f sec\n',ATtime./timing_trials);
fprintf('ABB  = %f sec\n',ABBtime./timing_trials);
fprintf('ATBB = %f sec\n',ATBBtime./timing_trials);
fprintf('-------------------\n');



