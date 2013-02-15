% fft1d_f: implementation of dense FFT-based sampling operator(e.g., partial FFT, scrambled FFT etc.); 
% also refer to fft1d_t: implementation of the transpose of dense
% FFT-based sampling operator; 
% Written by Thong Do, JHU, 2007
% Revised by Lu Gan, April. 2008

function b = fft1d_f(x, OMEGA, p)
% Input parameters:
    % x: Original signal vector;
    % OMEGA: a vector with K/2 frequency indices choosen uniform at random;
    % p: pre-randomized operator;
        %p is either a random permutation of [1:N] (for scrambling);
        %or a Bernoulli vector with 1, -1 entries (for random sign flipping); 
% Output Parameter: 
    % b: Sampled vector with length of K;

    % 1d DFT
    N = length(x);
    if max(p)>1,
        % pre-randomizing using permutation vector
        fx = 1/sqrt(N)*fft(x(p));
    elseif max(p) == 1
        % pre-randomizing using Bernoulli vector
        x = x.*p;
        fx = 1/sqrt(N)*fft(x);
    end
  b = [sqrt(2)*real(fx(OMEGA)); sqrt(2)*imag(fx(OMEGA))];