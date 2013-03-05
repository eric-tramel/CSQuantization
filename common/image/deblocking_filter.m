function Xrec = deblocking_filter(X,block_dim,cut_radius,method)
%DEBLOCKING_FILTER Remove block artifacts from an image
%
% function Xrec = deblocking_filter(X,block_dim,cut_radius)
%
% This function will remove blocking artifacts from an image 
% assuming the block dimensions given. Horizontal and 
% vertical strips of width 2*cut_radius will be replaced with
% an inpainting method. 
%
%	Inputs:
%
%	Outputs:
%
%
% Written by: Eric W. Tramel, Ph.D.
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

if nargin < 4
	method = 0;
end

% Image size
imsize = size(X);

% Set 'block' width
block_rows = imsize(1)./block_dim(1);
block_cols = imsize(2)./block_dim(2);


% Cut out block edges
for i = block_dim(2):block_dim(2):(imsize(2)-1)
	X(1:imsize(1),(i-cut_radius+1):(i+cut_radius)) = NaN;
end

for i = block_dim(1):block_dim(1):(imsize(1)-1)
	X((i-cut_radius+1):(i+cut_radius),1:imsize(2)) = NaN;
end

% Inpainting
Xrec = inpaint_nans(X,method);

