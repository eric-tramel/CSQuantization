function test2d_2(imgfile)
% This test compares the rate-distortion curve for the following
% sampling operators:
% 
%  1. Scrambled FFT (a dense operator), 
%  2. scrambled 32x32 block Hadamard transform (a highly sparse sampling
%  operator); 
%  3. Block DCT with random sign reversal of input signals (a sparse operator with full streaming capability); 

% Inputs: 
%   imgfile: input filename;
%   Example: test2('boat256.bmp');


img=imread(imgfile);
[row_num,col_num]=size(img);
alpha=0.05:0.1:0.75;
K=alpha*row_num*col_num;
N=row_num*col_num;
% Scrambled FFT;
[rec_sfft,psnr_sfft]=fast_cs2d(imgfile, K, 'FFT', 0);
% Scrambled Block Walsh-Hadamard Transform;
[rec_sbhe,psnr_sbhe]=fast_cs2d(imgfile, K, 'BWHT', 0,32);
% Block DCT + random sign flipping;
[rec_bdct,psnr_bdct]=fast_cs2d(imgfile, K, 'BDCT', 1,row_num);

close all;
plot(K/N,psnr_sfft,'b-',K/N,psnr_sbhe,'r--',K/N,psnr_bdct,'g--');
legend('Scrambled FFT','Scrambled BWHT','BDCT+Sign flipping');
xlabel('Sampling Rate');
ylabel('PSNR (dB)');
message=sprintf('R-D curve of %s',imgfile);
title(message);
