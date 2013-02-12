function run_cspl_ddwt(filename)

addpath(genpath('C:\Dropbox\work\Common\Sequences'));
addpath(genpath('C:\Dropbox\work\Common\l1magic'));
addpath(genpath('C:\Dropbox\work\Common\BCS-SPL'));
addpath(genpath('C:\Dropbox\work\Common\WaveletSoftware'));

%%% for test
filename = 'cardiac.000';
% filename = 'brain_sagittal';
%  subrates= 0.1:0.1:0.5;
subrates = 0.3;
num_trials = 1;
%%% for test

nonNegative = 1;
num_levels = 6;
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
      
      tic
      reconstructed_image = CS_PL_DDWT_Decoder(y, Phi, Phi_t, num_rows, ...
        num_cols, nonNegative, num_levels);
      toc
      
      reconstructed_image = reconstructed_image + im_mean;
      psnrs(j,i) = PSNR(original_image, reconstructed_image)
    end
 end
mean_psnrs = mean(psnrs,1);
% imwrite(uint8(reconstructed_image),[filename '_csspl_ddwt_radial.' num2str(subrates) '.pgm'])
% filename_save = [ './results/' filename '_csspl_ddwt.mat'];
% save(filename_save, 'psnrs','mean_psnrs');
% figure;imagesc(reconstructed_image);colormap gray;