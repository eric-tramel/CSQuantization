function test2d_1(imgfile,K)
% This example provides a quick comparison of reconstruced images offered
% by the following measurement operators: 
%  1. Scrambled FFT (a dense operator), 
%  2. scrambled 32x32 block Hadamard transform (a highly sparse sampling
%  operator); 
%  3. Block DCT with random sign reversal of input signals (a sparse operator with full streaming capability); 

% Input: 
%   imgfile: input filename;
%   K: No. of Samples; 

% Example: test2d_1('lena256.bmp',25000);
tic;
img=imread(imgfile);
% get the row and column dimension of the input image;
[r,c]=size(img);

% Scrambled FFT:
[rec_sfft,psnr_sfft]=fast_cs2d(imgfile, K, 'FFT', 0);

% Scrambled Block Walsh-Hadamard Transform;
[rec_sbhe,psnr_sbhe]=fast_cs2d(imgfile, K, 'BWHT', 0,32);

% Block DCT with random sign reversal;
% Block size=row number of input image;
[rec_bdct,psnr_bdct]=fast_cs2d(imgfile, K, 'BDCT', 1,r);
toc;
close all;
% Show the reconstructed image and the PSNR values:
figure(1);
imshow(rec_sfft);
message=sprintf('Scrambled FFT, PSNR=%6.2f dB',psnr_sfft);
title(message);

figure(2);
imshow(rec_sbhe);
message=sprintf('Scrambled 32x32 Block WHT, PSNR=%6.2f dB',psnr_sbhe);
title(message);

figure(3);
imshow(rec_bdct);
message=sprintf('Block DCT+ random sign flipping, PSNR=%6.2f dB',psnr_bdct);
title(message);


