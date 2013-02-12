%
% function s = SAM(x, x_hat)
%
%   This function calculates the angle, in degrees, between
%   vectors x and x_hat. SAM denotes "spectral angle mapper."
%

%
% CPPCA: Compressive-Projection Principal Component Analysis
% Copyright (C) 2007-2009  James E. Fowler
% http://www.ece.mstate.edu/~fowler
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%


function s = SAM(x, x_hat)

if (~isvector(x) || ~isvector(x_hat))
  error('x and x_hat must be vectors');
end

x = x(:);
x_hat = x_hat(:);

if (norm(x) == 0)
  s = 0;
else
  if (norm(x_hat) == 0)
    s = 0;
  else
    s = real(acos(x' * x_hat / norm(x) / norm(x_hat))) / pi * 180;
  end
end
