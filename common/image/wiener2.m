function f = wiener2(x, nhood)
% Wiener -- Wiener filter for 2D images variance is calculated using local
% mean and variance. Symetric extention is used for the boundary.
% Usage
%   f = Wiener2(x,nhood)
% Inputs
%   x           input image
%   nhood		neighborhood size
% Outputs
%   f		    denoised image
%
if (nhood(1) ~= 3) || (nhood(2) ~= 3)
  error('only support [3 3] for now');
end
if norm(x) == 0;
  f = zeros(size(x));
else
  % symetric extension 1 pixel around the boundary
  x = symextend(x,1);
  
  localMean = convn(ones(nhood), x) / prod(nhood);
  localVar = convn(ones(nhood), x.^2) / prod(nhood) - localMean.^2;
  
  localMean = localMean(2:end-1,2:end-1);
  localVar = localVar(2:end-1,2:end-1);
  
  noise = sum(x(:), [], 'double') / numel(x);
  
  f = x - localMean;
  x = localVar - noise;
  x = max(x, 0);
  localVar = max(localVar, noise);
  f = f ./ localVar;
  f = f .* x;
  f = f + localMean;
  
  % extract original data
  f = f(2:end-1,2:end-1);
end
function y = symextend(x,Nnum)
y = [fliplr(x(:,1:Nnum)) x x(:,end:-1:end-Nnum+1)];
y = [flipud(y(1:Nnum,:)); y ;y(end:-1:end-Nnum+1,:)];

