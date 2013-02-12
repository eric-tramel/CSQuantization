% 
% function [y Phi] = MC_BCS_SPL_Encoder(frames, subrate, key_subrate)
% 
%   This function performs BCS projections of each block of
%   the group of frames. Based on assumption that the key frame subrate
%   is larger than the non-key frame subrate, the program returns Phi 
%   whose rows are determined by key frame. The number of columns of 
%   the projection matrix, Phi, is determined by the size of the blocks
%   into which current_image is partitioned. Actual Projection is 
%   performed by BCS_SPL_Encoder frame by frame. The projections are 
%   returned as the columns of y.
%
%   See:
%     S. Mun and J. E. Fowler, "Residual Reconstruction For Block-
%     Based Compressed Sensing of Video," to appear in the Data
%     Compression Conference, 2011
%
%   Originally written by SungKwang Mun, Mississippi State University
%

%
% MC-BCS-SPL: Motion Compensated Block Compressed Sensing 
%               - Smooth Projected Landweber
% Copyright (C) 2011  James E. Fowler
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

function [y Phis] = MC_CS_SPL_Encoder(frames, key_subrate, subrate)

[num_rows num_cols] = size(frames{1});
num_frames = length(frames);

% key frame
L = round(165*key_subrate^2 + 218*key_subrate + 3);
[MM,Mh,mh,OMEGA] = LineMask(L,num_rows);

Phi_key = @(z) A_fhp(z, OMEGA);
Phi_key_t = @(z) At_fhp(z, OMEGA, num_cols);

y = cell(1,num_frames);

y{1} = Phi_key(frames{1}(:));

% non-key
L = round(165*subrate^2 + 218*subrate + 3);
[MM,Mh,mh,OMEGA] = LineMask(L,num_rows);
Phi = @(z) A_fhp(z, OMEGA);
Phi_t = @(z) At_fhp(z, OMEGA, num_cols);

for i = 2 : num_frames
    y{i} = Phi(frames{i}(:));
end

y{end} = Phi_key(frames{end}(:));

Phis = {Phi_key Phi_key_t Phi Phi_t key_subrate subrate};
