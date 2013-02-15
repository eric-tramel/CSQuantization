% blk_f1d: implementation of the random block-based sensing operator;
% refer to blk_t1d: implemenatation of the transpose of random block-based sensing
% operator;
% Written by Thong Do, JHU, 2007,Feb. 2008
% Revised by Lu Gan, Brunel University, April. 2008

function b = blk_f1d(x, select_vect, rand_vect,trans_type,blk_size)
% Input parameters
% x: Original signal vector;
% select_vect: vector that contains measurement indices; 
% rand_vect: vector that randomizes samples 
    %it is either a random permutation of [1:N] 
    %or a Bernoulli vector with 1, -1 entries;  
% trans_type: block operator, either 'DCT' or 'WHT';
% blk_size: block size of the input signal; 

% Output:
% b: Sampled signal vector;

% Examples: 
% 1. Scrambled Block Hadamard Ensemble with block size of 32;
%   b = blk_f1d(x, select_vect, randperm(1:length(x)),'BWHT',32);
% 2. Block DCT-based measurement with random flipping of signal signs (block
% size=512); 
%   b = blk_f1d(x, select_vect, randperm(1:length(x)),'BDCT',512);


% Find out the length of the input signal;
N = length(x);

% Check the block size;
if mod(N,blk_size)>0
    error('length(x)/blk_size must be an integer');
end

% Check the trans_type;
trans_type=upper(trans_type);

if strcmp(trans_type,'BDCT')
    Phi_B=dct(eye(blk_size));
elseif strcmp(trans_type,'BWHT')
    Phi_B=hadamard(blk_size)/sqrt(blk_size);
else
     error('Block trans_type should be either DCT or WHT');
end


% Pre-randomize the signal;
if max(rand_vect)>1
    % Random permutation; 
    x(rand_vect) = x;
elseif max(rand_vect)==1
    % Random sign flipping; 
    x = x.*rand_vect;
end

% Reshape x into a 2D matrix for paralell processing; 
x_blk = reshape(x,blk_size,floor(N/blk_size));
% Transform of the signal; 
fx = Phi_B*x_blk;
% Vectorize the signal;
fx = fx(:);
% Select the samples; 
b = fx(select_vect);