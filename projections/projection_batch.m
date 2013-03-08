function Y = projection_batch(A,X,B,N)
%PROJECTION_BATCH Applies projection A to the columns of X
%
% function A = projection_batch(A,x,B,M)
%
%  Applies projection A to the columns of X, assuming that A
% is a function handle. If A is not a function handle, then this
% function will break. 
%
%	Inputs:
%		A 		Function handle for projection.
%		X 		A NxB matrix of N dimensional vectors to apply A to.
%  		B,N 	Dimensions of X (optional)	
%
%	Outputs:
%		Y 		A ?xB matrix of projections.
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


% TODO: There should be a faster way to do this batch projection than
% through a loop. Need to look into methods for applying a 
% function along columns.
%
% It seems that the fastest approach may be to build a batch version of the
% blk_t1d, blk_f1d functions to use with block based projections. The
% reasoning for this is that the majority of time is taken up by the
% 'hadamard' functions which need to be regenerated on each blk_f1d/blk_t1d
% call. In batch mode these would only be generated once.


%% Looping Approach
if nargin > 2
	Y = zeros(N,B);
end

for i=1:B
	Y(:,i) = A(X(:,i));
end

%% Cell Approach
% Found in practice that this approach is not any faster than looping
% X = num2cell(X,1);
% Y = cell2mat(cellfun(A,X,'UniformOutput',false));






