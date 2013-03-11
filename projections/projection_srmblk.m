function [A AT B BT] = projection_srmblk(subrate,N,trans_mode,blksize,is_blocked,imsize,block_dim)
%PROJECTION_SRMBLK generate funciton handles for SRM-block random projection.
%
% function [A AT B BT] = projection_srmblk(M,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb)
%
% Creates a set of funciton handles for both forward and backward random projections using 
% structured random matrices. This function uses the block mode SRMs. Utilizes the SRM
% toolbox software written by Do, Gan, & Tran. 
%
%	Inputs:
%
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

%% Input Handling
% Check to see if we are in image or vector projection mode
Nb = 0;
if nargin < 5
	is_blocked = 0;
	imsize = [];
	block_dim = [];
end

%% Block-based settings
if is_blocked
    N = block_dim(1)*block_dim(2);
    blocks = imsize ./ block_dim;
    Nb = blocks(1)*blocks(2);
end

%% Get Measurements
M = round(subrate*N);

if M > N
    rescale = M ./ N;
else
    rescale = 1;
end

%% Forward Projection
% Split up projection sets to allow for Subrates > 1
Mleft = M;
i = 1;
while Mleft >= N
	% Each 'sub' projection has a dimensionality one less than the original 
	% signal dimensionality. This is to avoid any odd processing errors from 
	% the SRM code for the (Subrate = 1.0) case.
	thisM = (N - 1);
	% Generate sub projection
	[B{i} BT{i}] = single_projector_set(thisM,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb);
	i = i + 1;
	Mleft = Mleft - thisM;
end

% If there are left over measurements
if Mleft ~= 0
	[B{i} BT{i}] = single_projector_set(Mleft,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb);
end

% Set up forward transform
% To speed up this process, we will perform a check to see if we indeed
% have more than one projection
if length(B) > 1
    A = @(z) srm_forward_batch(B,z);
else
    A = @(z) B{1}(z);
end


%% Backward Projection
% Set up transpose transform. Needs to account for whether or not the 
% measurements have to be rearranged in block-based mode.
% To do this, we need to sum together all the individual projections
if is_blocked
    AT = @(z) srm_transpose_batch(BT,z,N*Nb,Nb,i,rescale);
else
    AT = @(z) srm_transpose_batch(BT,z,N,1,i,rescale);
end

%% Helper functions
%----------------------------------------------------
function [A AT] = single_projector_set(M,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb)
% This is the function which actually generates the random projection function handles.
rand_vect = randperm(N)';
select_vect = randperm(N);
select_vect = select_vect(1:M);
Phi   = @(z) blk_f1d(z,select_vect,rand_vect,trans_mode,blksize);
PhiT  = @(z) blk_t1d(z,N,select_vect,rand_vect,trans_mode,blksize);

if is_blocked
    A = @(z) srm_forward(Phi,z,M,imsize,block_dim,Nb);
    AT = @(z) srm_transpose(PhiT,z,M,N,block_dim,imsize,Nb);
else
	A = Phi;
	AT = PhiT;
end

function y = srm_forward(Phi,x,M,imsize,block_dim,Nb)
% Using a separate function to avoid in-line definition mess.
y = reshape(x,imsize);
y = im2col(y,block_dim,'distinct');
y = projection_batch(Phi,y,Nb,M);
y = y(:);

function x = srm_transpose(PhiT,y,M,N,block_dim,imsize,Nb)
x = reshape(y,[M,Nb]);
x = projection_batch(PhiT,x,Nb,N);
x = col2im(x,block_dim,imsize,'distinct');
x = x(:);

% function y = batch_projection(A,x,M,B)
% % There should be a faster way to do this batch projection than
% % through a loop. Need to look into methods for applying a 
% % function along columns;
% y = zeros(M,B);
% for i=1:B
% 	y(:,i) = A(x(:,i));
% end

function x = srm_transpose_batch(AT,y,N,Nb,num_projs,rescale)
x = zeros(N,1);
i = 1;
while i < num_projs
	x = x + AT{i}(y(1:(N-Nb)));
	y = y((N-Nb+1):end);
	i = i + 1;
end
x = x + AT{i}(y);
x = x ./ rescale;

function y = srm_forward_batch(A,x)
% Assuming A is a cell array
y = [];
for i=1:length(A)
    y = vertcat(y,A{i}(x));
end