% This function implements fast spatial-domain compressive sampling of 2D
% images through structurally random matrices (SRM). Details of SRM can be
% found in the following paper:  

% T. Do, T. D. Tran and L. Gan "Fast compressive sampling with structurally random matrices", ICASSP'08

% The l_1 reconstructed method is based on the GPSR_BB package at www.lx.it.pt/~mtf/GPSR/;

function [rec_img,psnr_val]=fast_cs2d(im, K, trans_mode, rand_type, blk_size);
%  ===== Required inputs ============= 
% im : filename of the input image

% K: number of measurements
  % K can be either a vector or a scalar; 
  % If K is a scalar, it will show the reconstructed image;
  % If K is a vector, the function will produce a R-D curve;  

% trans_mode: measurement matrix type
    % 'PFFT': Partial FFT without pre-randomization;
    % 'FFT': FFT with pre-randomization;
    % 'BDCT': Block DCT with pre-randomization;
    % 'BWHT': Block Walsh-Hadamard Transform (WHT) with pre-randomization;

% rand_type: type of pre-randomization operator. 
    %0=random permutation (global model);
    %1=random flipping the signs of the signal(local model); 

%  ===== Optional Input =============
% blk_size: block size of the measurement operator if the trans_mode is either 
%           'BDCT' or 'BWHT'; 
%           Default: 32 for random permutation (rand_type=0) and row_num of input image
%           for random flipping the signs (rand_type=1);

% **********************************************************************

% ===== Outputs =============
% rec_img: reconstructed image if K is a scalar; []if K is a vector; 
% psnr_val: the PSNR of reconstruced image;  

% For Examples on how to use this program, please refer to test2d_1.m,
% test2d_2.m and test2d_3.m;


% Test the input parameters:
if nargin<4,
    help fast_cs2d;
    error('Wrong Number of Input Parameters');
end


% Read the input image;            
x = double(imread(im));
% keep an original copy of the input signal
x0 = x;


[m n] = size(x);
% Total number of pixels;
N = m*n;
x = x(:);
% Substract the mean of the input signal; 
xmean = mean(x);
x = x-xmean;

% Set the default of blk_size for BDCT or BWHT if it is not specified; 
if (strcmp(trans_mode,'BDCT') | strcmp(trans_mode,'BWHT'))& (nargin<5)
   if rand_type==0,
       blk_size=32;
   else
       blk_size=m;
   end
end

% Sparsifying transform: 9-7 Wavelet 
[h0,h1,f0,f1] = filter9_7();
L = floor(log2(m))-3;               % Level of decomposition

% Initialize the output parameters;
rec_img=[];
psnr_val=[];

% Choose an arbitrary random seed; 
% User can change it
Perm_state = 3587642;
rand('state', Perm_state);
  
% Define the random vector for the input samples: 

% Random permutation;
if rand_type == 0
    if strcmp(trans_mode,'PFFT'),
        rand_vect=1:N; %No random permutation for PFFT;
    else
   % other modes: randomize samples use permuation vector
        rand_vect = randperm(N)';
    end
% Random filipping the sign;    
elseif rand_type == 1
    if strcmp(trans_mode,'PFFT'),
        rand_vect=ones(size(x)); %No sign flipping for PFFT;
    else
        % Other modes: randomize samples using a Bernoulli vector
        rand_vect = 2*round(rand(N,1))-1;
    end
end % of if rand_type

% Main loop to test the reconstruction performance for each K(i);
 for i =1:length(K),
       Ki=round(K(i)); 
       
       % Dense FFT-based Operators;
        if strcmp(trans_mode,'PFFT')| strcmp(trans_mode,'FFT')
            % Define selected samples
            select_vect = randperm(round(N/2)-1)+1;
            select_vect = select_vect(1:round(Ki/2))';
            % Define Sampling Operator;
            Phi = @(z) fft1d_f(z, select_vect, rand_vect);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) fft1d_t(z, N, select_vect, rand_vect);  
            
        % Block-Based Operators:
        elseif strcmp(trans_mode,'BDCT') | strcmp(trans_mode,'BWHT')
                        
            % Define selected samples
            select_vect=randperm(N);
            select_vect = select_vect(1:Ki)'; 
            
            % Define Sampling Operator;    
            Phi = @(z) blk_f1d(z, select_vect, rand_vect, trans_mode, blk_size);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) blk_t1d(z, N, select_vect, rand_vect,trans_mode, blk_size); 
            
        end

   % Define the equivalent sampling operator for the wavelet coefficients; 
   B = @(z) A_idwt1d(Phi,z,f0,f1,L,m,n);
   Bt = @(z) At_dwt1d(Phi_T,z,h0,h1,L,m,n);
   
    % getting measurements
     y = Phi(x);
   
    % Reconstruction. Using GPSR_BB modules
    tau = norm(y,'fro')/sqrt(Ki)/16;   
    [alp,alp_debias,objective,times,debias_start,mses]= ...
         GPSR_BB(y,B,tau,...
         'AT', Bt,'Debias',1,'Initialization',0,...
         'StopCriterion',1,'ToleranceA',0.001,'ToleranceD',0.00005);
     
   % Transform from the WT domain to the spatial domain  
    alp_debias = reshape(alp_debias,m,n);
    xr = idwt2d(alp_debias,f0,f1,L);
    % add mean to the reconstruction
    xr = xr+xmean;
    psnr_val  = [psnr_val,psnr(x0,xr)];
 end
 
% Retrun the reconstructed image if K is a scalar; 
if length(K)==1,
    rec_img=uint8(xr);
end


     