% 
% function [reconstructed_frames reconstructed_frames_indep] = ...
%     MC_BCS_SPL_Decoder_ForwardBackward(y, Phi, num_rows, num_cols)
% 
%   This function performs consecutive reconstruction which starts from
%   the key frame of GOP which normally sampled at high subrate. 
%   In this way, reconstructed current frame will be reference frame for
%   the next frame. This function returns reconstructed frames by MC-BCS-SPL
%   and independent reconstruction by BCS-SPL. 
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

function [reconstructed_frames recon_3ds] = ...
    MC_CS_PL_Decoder_ForwardBackward(y, Phis, num_rows, num_cols, recon_3d, recon_3ds, frames)

  num_frames = length(y);
  reconstructed_frames = cell(1,num_frames);
  reconstructed_frames{1} = recon_3ds{1};
  
  for i = 2:num_frames
    recon_3ds_without_current = {recon_3ds{1:i-1} recon_3ds{i+1:end}};
    reconstructed_frames{i} = ...
      MC_CS_PL_Decoder(reconstructed_frames{i-1}, recon_3d{i-1}, y{i}, Phis, num_rows, num_cols, recon_3ds_without_current);
    recon_3ds = {recon_3ds{1:i-1} reconstructed_frames{i} recon_3ds{i+1:end}};
  end
