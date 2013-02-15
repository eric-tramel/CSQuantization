function test2d_3(img_file)
% This example aims to re-produce Figure 1 and Figure 2 in the following
% paper: 
% T. Do, T. D. Tran and L. Gan "Fast compressive sampling with structurally
% random matrices", ICASSP'08

% To reproduce Figure 1: test2d_3('lena512.bmp');
% For figure 2:test2d_3('boat512.bmp');
tic;
samp_rate=0.05:0.1:0.75;
K=round(samp_rate*512*512);

% Benchmark: Wavelet domain measurement through partial FFT 
[rec_img,psnr_bf]=pfft_cswt2d(img_file, K);
N=512*512;

% Spatial domain measurements:
[rec_img,psnr_pfft]=fast_cs2d(img_file, K, 'PFFT', 0);
[rec_img,psnr_sfft]=fast_cs2d(img_file, K, 'FFT', 0);
[rec_img,psnr_dct512]=fast_cs2d(img_file, K, 'BWHT', 1,512);
[rec_img,psnr_dct32]=fast_cs2d(img_file, K, 'BWHT', 0,32);
[rec_img,psnr_wht512]=fast_cs2d(img_file, K, 'BDCT', 1,512);
[rec_img,psnr_wht32]=fast_cs2d(img_file, K, 'BDCT', 0,32);

close all;
figure;
plot(K/N,psnr_bf,'r-',K/N,psnr_pfft,'g--',K/N,psnr_sfft,'b^',K/N,psnr_dct512,'bo',K/N,psnr_dct32,'ms',K/N,psnr_wht512,'b*',K/N,psnr_wht32,'m+');
legend('BF','PFFT','SFFT','DCT512','DCT32','Ha512','Ha32');
title('R-D performance');
xlabel('Sampling Rate');
ylabel('PSNR (dB)');
toc;