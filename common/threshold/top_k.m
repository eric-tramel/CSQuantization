function x = top_k(x,k)
% x = top_k(x,k)
% Thresholding function which zeros all values in x except the top K with
% largest magnitude.

[toss,idx] = sort(abs(x),'descend');
x(idx(k+1:end)) = 0;