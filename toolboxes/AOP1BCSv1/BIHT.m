function [x iter] = BIHT(b, A, x, K, maxiter, alpha) 
%%
% BIHT is the function used to recover signals from 1-bit measurements in 
% 1-bit compressive sensing framework. It uses BIHT algorithm proposed in
% the reference.
%
%                   Min.  phi(b,Ax)
%                   s.t.  |x|_0 <= K,  |x|_2 = 1,
%
%
% Inputs are 
%            b:   The 1-bit measurements (1 or -1), maybe corrupted                       
%            A:   The measurement system (M-by-N matrix)
%            x:   Initial guess of the signal
%            K:   sparsity of the signal (only K components of x are
%                 non-zero)
%      maxiter:   The maximum total number of iterations
%        alpha:   A parameter used in the gradient descent step                       
%
% Outputs are
%            x:   The recovered signal
%         iter:   Total number of iterations used
%         
%
% Author:    Ming Yan (basca.yan@gmail.com) and Yi Yang (yyang@math.ucla.edu.cn)
%   Date:    2011-10-24(UCLA)
% 
% Reference: 
%          Laurent Jacquesy, Jason N. Laskaz, Petros T. Boufounosx, and
%          Richard G. Baraniukz, Robust 1-Bit Compressive Sensing via
%          Binary Stable Embeddings of Sparse Vectors
%
%%

    htol    = 0;
    hd      = Inf;
    iter    = 0;

    while(htol < hd) && (iter < maxiter)        
        x   = x + alpha * A' * (b - sign(A*x));
        [~, index] = sort(abs(x),'descend');
        x(index(K+1:end)) = 0;
        hd   = nnz(b - sign(A*x));
        iter = iter + 1;
    end
end