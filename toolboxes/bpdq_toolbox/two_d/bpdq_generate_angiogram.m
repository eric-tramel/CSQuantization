function result = bpdq_generate_angiogram(m,n,num,seed)
% bpdq_generate_angiogram : Generate artificial Angiogram signal
%
% function result = bpdq_generate_angiogram(m,n,num,seed)
%
% Inputs:
% m,n - size of desired synthetic angiogram
% num - number of ellipses
% seed - random seed
% Outputs:
% result - computed synthetic angiogram
%
% Adapted from code from Ewout van den Berg and Michael P. Friedlander
%
% GENERATEANGIOGRAM  Generate artificial Angiogram signal.
%
%   GENERATEANGIOGRAM(M,N) generates an M by N image filled with 18
%   non-overlapping ellipses.
%
%   Based on the SparseMRI implementation by M. Lustig, 2007.
%
%   References
%      M. Lustig, D.L. Donoho, and J.M. Pauly, Sparse MRI: The
%         application of compressed sensing for rapid MR imaging,
%         Submitted to Magnetic Resonance in Medicine, 2007.
%
% This file is part of Sparco Toolbox (http://www.cs.ubc.ca/labs/scl/sparco/)
% Copyright (C) 2008, Ewout van den Berg and Michael P. Friedlander   
% (See the notice at the end of the file.)


% Initialise the random number generator
  if~exist('seed','var')
    seed=16000;
  end
  
  rand('twister',seed);

% Settings
imgSize    = [m,n];
nMagnitude = 3;
%sizes      = [1,3,6];
sizes      = [1,3,6];
magnitudes = (1:nMagnitude) / nMagnitude;

% Initialise the vectors
n  = round(length(sizes) * (length(sizes) + 1) / 2);
sx = []; sy = [];
mg = zeros(n * nMagnitude,1);

% Set up the x- and y-radii and intensity values
mg = repmat(magnitudes,n,1);
for i=1:length(sizes)
  sx = [sx; repmat(sizes(i),   i,nMagnitude)];
  sy = [sy; repmat(sizes(1:i)',1,nMagnitude)];
end
sx = sx(:) * (imgSize(2) / 100.0);
sy = sy(:) * (imgSize(1) / 100.0);
mg = mg(:);

sx=sx(1:num);
sy=sy(1:num);

% ---------------------------------------------------------------------
% Generate set of non-overlapping ellipsoids.
% Allow at most 'maxfail' occurences of overlapping ellipsoids
% to avoid potential infinite looping in dense settings.
% ---------------------------------------------------------------------

maxfail = 1000;
m       = imgSize(1);
n       = imgSize(2);
img     = zeros(m,n); % Image
ell     = zeros(m,n); % Single ellipsoid
[x,y]   = meshgrid((0:n-1)-round(n/2), (0:m-1)-round(m/2));

for i=1:length(sx)
   angle = 2*pi*(rand(1)-0.5);
   while (maxfail > 0)
      ell = ellipsoid_intrnl(sx(i),sy(i),angle);
      if sum(any(ell .* img)) == 0, break; end; % No overlap
      maxfail = maxfail - 1;
   end
   img = img + mg(i) * ell;
end

result = img;


function res = ellipsoid_intrnl(xRadius,yRadius,angle)

   b  = sqrt(2 / (1 + xRadius^2 / yRadius^2));
   a  = b * xRadius / yRadius;
   xr = a * (x * cos(angle) - y * sin(angle));
   yr = b * (y * cos(angle) + x * sin(angle));
   mx = max(abs(xr(:)));
   my = max(abs(yr(:)));

   res = 0;
   while (sum(res(:)) == 0) % While empty plot
      xc  = 1.1 * (rand(1,1)-0.5) * mx;
      yc  = 1.1 * (rand(1,1)-0.5) * my;
      res = (((xr-xc).^2 + (yr-yc).^2) < yRadius^2);
   end
end

end % program

% The Sparco Toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%  
% The Sparco Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License
% along with Sparco Toolbox. If not, see <http://www.gnu.org/licenses/>.
