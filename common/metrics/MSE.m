function mse = MSE(x,y)
% mse = MSE(x,y)
%
% Calcualte the root meat square between two vectors. 
x = x(:); y = y(:);
N = length(x);

if N ~= length(y)
    error('rms:DimMismatch','Input dimensionalities must match.');
end

r = x - y;
mse = r'*r ./ N;
