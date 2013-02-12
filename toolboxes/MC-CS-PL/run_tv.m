function run_csspl_ddwt(filename)

% for medical images, mean value should be subtracted.
% PSNR difference is about 5dB

addpath(genpath('C:\Dropbox\work\Common\Sequences'));
addpath(genpath('C:\Dropbox\work\Common\SRM'));
addpath(genpath('C:\Dropbox\work\Common\l1magic'));
addpath(genpath('C:\Dropbox\work\Common\BCS-SPL-DPCM'));
addpath(genpath('C:\Dropbox\work\Common\BCS-SPL'));
addpath(genpath('C:\Dropbox\work\Common\WaveletSoftware'));

addpath(genpath('/home/sm655/Common/Sequences'));
addpath(genpath('/home/sm655/Common/BCS-SPL'));
addpath(genpath('/home/sm655/Common/MS-BCS-SPL'));
addpath(genpath('/home/sm655/Common/WaveletSoftware'));
addpath(genpath('/home/sm655/Common/BCS-SPL-DPCM'));
addpath(genpath('/home/sm655/Common/l1magic'));

%%% for test
% filename = 'cardiac.000';
filename = 'brain_sagittal';
%  subrates= 0.1:0.1:0.5;
subrates = 0.3
num_trials = 1;
%%% for test

global original_image
if strcmp(filename,'phantom')
  % Phantom
  n = 256;
  original_image = phantom(n);
else
  original_filename = [filename '.pgm'];
  original_image = double(imread(original_filename));
end
[num_rows num_cols] = size(original_image);

% N = num_rows*num_cols;

% im_mean = mean(original_image(:));
im_mean = 0;
x = original_image - im_mean;

for i = 1:length(subrates);
    subrate = subrates(i);

    L = round(165*subrate^2 + 218*subrate + 3)

    for j = 1:num_trials
      % Fourier samples we are given
      [MM,Mh,mh,mhi] = LineMask(L,num_rows);
      OMEGA = mhi;
      
      Phi = @(z) A_fhp(z, OMEGA);
      Phi_t = @(z) At_fhp(z, OMEGA, num_cols);

%       [pdf,val] = genPDF([num_rows num_cols],power,subrate,2,0,0);
%       [mask1,stat,NN] = genSampling(pdf,100,10);
%       mask1=fftshift(mask1);

      y=Phi(x(:));
      

 % min l2 reconstruction (backprojection)
      xbp = Phi_t(y);

      % recovery
      tic
      xp = tveq_logbarrier(xbp, Phi, Phi_t, y, 1e-3, 10, 1e-8, 200);
      reconstructed_image = reshape(xp, num_rows, num_cols);
      toc
      
      
      reconstructed_image = reconstructed_image + im_mean;
      psnrs(j,i) = PSNR(original_image, reconstructed_image)
    end
 end
mean_psnrs = mean(psnrs,1);
imwrite(uint8(reconstructed_image),[filename '_tv_radial.' num2str(subrates) '.pgm'])
% filename_save = [ './results/' filename '_csspl_ddwt.mat'];
% save(filename_save, 'psnrs','mean_psnrs');
% figure;imagesc(reconstructed_image);colormap gray;