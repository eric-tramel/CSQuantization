function [x iter Loc] = BIHT2_AOP_flip(b, A, x, K, L, miter, maxiter, alpha) 
%%
% BIHT2_AOP_flip is the function used to recover signals from 1-bit measurements in 
% 1-bit compressive sensing framework. It detects the locations of sign flips 
% and defines new measurements by flipping the signs of the input measurements 
% at these locations. In each step it utilizes BIHT-L2 algorithm with this 
% new measurements and the measurement system.
%
%                   Min.  Lambda.*phi(b,Ax)
%                   s.t.  |x|_0 <= K,  |x|_2 = 1,
%                         sum(1-Lambda_i) <= L, Lambda_i = 0 or 1.  
%
%
% Inputs are 
%            b:   The 1-bit measurements (1 or -1), with L measurements having sign flips                       
%            A:   The measurement system (M-by-N matrix)
%            x:   Initial guess of the signal
%            K:   sparsity of the signal (only K components of x are
%                 non-zero)
%            L:   The maximum number of sign flips in the measurements
%        miter:   The number of inner iterations (usually 1 is enough)
%      maxiter:   The maximum total number of iterations
%        alpha:   A parameter used in the gradient descent step                       
%
% Outputs are
%            x:   The recovered signal
%         iter:   Total number of iterations used
%          Loc:   The position of mmeasurements detected as "correct"
%                                                                       
%
% Author:    Ming Yan (basca.yan@gmail.com) and Yi Yang (yyang@math.ucla.edu)
%   Date:    2012-3-26(UCLA)
%
% Reference:  M. Yan, Y. Yang and S. Osher, Robust 1-bit compressive
% sensing using adaptive outlier pursuit. 
%	      UCLA CAM report 11-71.
%
%%
    htol    = L;
    hd      = Inf;
    hd1     = Inf;
    iter    = 0;
    y_t     = A*x;
    % initial "correct" location
    Loc     = 1:size(A,1);
    % initial "corrected" measurements
    bn      = b;
    y       = bn.*max(bn.*y_t,0);
    
    while(htol < hd)&&(iter < maxiter)   
         % x-update with BITH-L2
        for inniter = 1:miter
            x       = x + alpha * A' * (y - A*x);    
            [~, index] = sort(abs(x), 'descend');
            x(index(K+1:end)) = 0;
            y_t     = A * x;
            y       = bn .* max(bn.*y_t, 0);
        end
        hd   = nnz(b - sign(y_t));
        % only update the Loc when fewer sign flips are found
        if hd < hd1
            hd1     = hd;
            % find the largest L elements of phi(b,Ax)
            [~, index] = sort(abs(y_t) .* max(-sign(b.*y_t), 0), 'descend');
            % update the "correct" location
            Loc     = 1:size(A,1);
            Loc(index(1:L)) = [];
            % correct the measurements by flipping the sign of b at these locations
            bn      = b;
            bn(index(1:L)) = -bn(index(1:L));
            y       = bn .* max(bn.*y_t, 0);
        end
        
        iter = iter + 1;
    end
end