% blk_t1d: implemenation of the transpose of the random block-based sensing operator;
% refer to blk_f1d: implemenatation of the random block-based sensing
% operator;
% Written by Thong Do, JHU, 2007
% modified by Lu Gan, Brunel University, April 2008

function x = blk_t1d(b, N, select_vect, rand_vect, trans_type, blk_size)
% Input parameters
% b: Sampled signal vector;
% N: length of the original (reconstructed) signal vector;
% select_vect: vector that contains measurement indices chosen uniform at
%              random;
% rand_vect: vector that randomizes samples
% rvec: rvec = 0 if rand_vect is permutation. rvec = 1 if rand_vect is Bernoulli
% Phi1: left multiplication matrix
% Phi2: right multiplication matrix
% Output
% x: (reconstructed) signal with length of N;

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


K = length(b);
fx = zeros(N,1);
fx(select_vect) = b(1:K);
% Reshape the vector into a 2D matrix for paralell processing;
fx = reshape(fx,blk_size,floor(N/blk_size));
x = Phi_B'*fx;
x = x(:);
if max(rand_vect) > 1
    % post-randomizing using permutation vector (conjugate of the pre-randomizing
    % operation)
    x = x(rand_vect);
elseif max(rand_vect) == 1
    % post-randomizing using Bernoulli vector (conjugate of the pre-randomizing
    % operation)
    x = x.*rand_vect;
end