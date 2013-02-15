% This function implements compressive sampling of wavelet coefficients of a 2D
% image through partial FFT measurement; It is presented here to serve as a benchmark for spatial 
% domain measurement operator;
% The l_1 reconstructed method is based on the GPSR_BB package;

function [rec_img,psnr_val]=pfft_cswt2d(im, K);
%  ===== Required inputs ============= 
% im : filename of the input image

% K: number of measurements
  % K can be either a vector or a scalar; 
  % If K is a scalar, it will show the reconstructed image;
  % If K is a vector, the function will produce a R-D curve;  


% **********************************************************************
% Outputs:
% rec_img: reconstructed image if K is a scalar; []if K is a vector; 
% psnr_val: the PSNR of reconstruced images;  

% Test the input parameters:
if nargin<2,
    help pfft_cswt2d;
    error('Wrong Number of Input Parameters');
end


% Read the input image;            
x = double(imread(im));
% keep an original copy of the input signal
x0 = x;


[m n] = size(x);
% Total number of pixels;
N = m*n;
% Substract the mean of the input signal; 
xmean = mean(x(:));
x = x-xmean;

% Sparsifying transform: 9-7 Wavelet 
[h0,h1,f0,f1] = filter9_7();
L = floor(log2(m))-3;               % Level of decomposition

% Get the wavelet transform coefficients:
xw = dwt2d(x,h0,h1,L);
xw = xw(:);

% Initialize the output parameters;
rec_img=[];
psnr_val=[];

% Main loop to test the reconstruction performance for each K(i);
 for i =1:length(K),
       Ki=round(K(i)); 
       
       % Define selected samples
         select_vect = randperm(round(N/2)-1)+1;
         select_vect = select_vect(1:round(Ki/2))';
         % Define Sampling Operator;
         Phi = @(z) fft1d_f(z, select_vect, 1:N);
         % Define the transpose of the Sampling Operator;
         Phi_T = @(z) fft1d_t(z, N, select_vect, 1:N);  
 
        % getting measurements
         y = Phi(xw);
   
    % Reconstruction. Using GPSR_BB modules
    tau = norm(y,'fro')/sqrt(Ki)/16;   
    [alp,alp_debias,objective,times,debias_start,mses]= ...
         GPSR_BB(y,Phi,tau,...
         'AT', Phi_T,'Debias',1,'Initialization',0,...
         'StopCriterion',1,'ToleranceA',0.001,'ToleranceD',0.00005);
     
   % Transform from the WT domain to the spatial domain  
    alp_debias = reshape(alp_debias,m,n);
    xr = idwt2d(alp_debias,f0,f1,L);
    % add mean to the reconstruction
    xr = xr+xmean;
    psnr_val  = [psnr_val,psnr(x0,xr)];
 end

figure; 

if length(K)==1,
    imshow(uint8(xr)); 
    message=sprintf('PFFT measurement in the wavelet domain; PSNR=%6.2f dB', psnr_val); 
    title(message);
    rec_img=uint8(xr);
else
      plot(K/N,psnr_val);
      ylabel('PSNR (dB)');
      xlabel('Rate K/N');
      message=sprintf('R-D curve of %s; PFFT measurement in the wavelet domain', im); 
      title(message);
end


     