function rms = RMS(x,y)
% rms = RMS(x,y)
%
% Calcualte the root meat square between two vectors. 
% N = length(x(:));
% 
% if N ~= length(y(:))
%     error('rms:DimMismatch','Input dimensionalities must match.');
% end
mse = MSE(x,y);
rms = sqrt(mse);