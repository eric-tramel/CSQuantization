% 
%function [reconstructed_center_frame recon_indep] = ...
%    MC_BCS_SPL_Decoder_Center(y, previous_frame, next_frame, Phi)
% 
%   This function performs MC-BCS-SPL method to reconstruct center frame
%   based on two reference frame in bidirecional way. Mulityhypothesis 
%   initialization includes independent reconstruction and two residual
%   reconstruction without motion. Two residual reconstruction of motion 
%   compensated frames composited by adjacent frames are averaged as final
%   reconstruction of center frame. 
%   This function returns reconstructed frame by MC-BCS-SPL
%   and independent reconstruction by BCS-SPL. 
%
%   See MC-BCS-SPL section on:
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

function reconstructed_center_frame = ...
    MC_CS_PL_Decoder_Center(y, previous_frame, next_frame, Phis, recon_3d, recon_3ds, center_frame)

subrate = Phis{6};
weight = subrate*2 - 0.2;
% weight = 0.5;
weight(weight>0.8) = 0.8;
weight(weight<0.2) = 0.2;

GOP_size = length(recon_3ds) - 1; 
ind_center = floor(GOP_size/2) + 1;

[num_rows num_cols] = size(previous_frame);

% Hypothesis Initialization
recon_res_zeroMV1 = ...
  MC_CS_PL_ResidualReconstruction(previous_frame, y, Phis, num_rows, num_cols);
recon_res_zeroMV2 = ...
  MC_CS_PL_ResidualReconstruction(next_frame, y, Phis, num_rows, num_cols);
recon_without_motion = (recon_res_zeroMV1 + recon_res_zeroMV2)/2;

reconstructed_current_frame = recon_3d*weight + recon_without_motion*(1-weight);

% Bidirectional iterative motion compensation

recon_3ds{ind_center} = reconstructed_current_frame;

num_iter = 1; 
accuracy = 'fp';
block_size = 32;
window_size = 7;
reference_frame = previous_frame;
recon_current_frame_forward = ...
    MC_CS_PL_IterativeReconstruction(reconstructed_current_frame, ...
    reference_frame, y, Phis, num_iter, block_size, window_size, recon_3ds, accuracy,weight); 

reference_frame = next_frame;
recon_current_frame_backward = ...
    MC_CS_PL_IterativeReconstruction(reconstructed_current_frame, ...
    reference_frame, y, Phis, num_iter, block_size,window_size, recon_3ds, accuracy,weight); 
% average result
reconstructed_center_frame = ...
    (recon_current_frame_forward + recon_current_frame_backward)/2;

% Auxiliary function
function reconstructed_frame = ...
  MC_CS_PL_IterativeReconstruction(reconstructed_frame, reference_frame, ...
  y, Phis, num_iter, block_size,window_size, recon_3ds, accuracy,weight)

%  global original_image
[num_rows num_cols] = size(reconstructed_frame);
for i = 1:num_iter
  motion_compensated_image1= memc_telescope(recon_3ds, reconstructed_frame, block_size, window_size);
  motion_compensated_image2 = motion_estimation_compensation(reference_frame, reconstructed_frame, ...
    block_size, window_size);
  %     motion_compensated_image2 = MC_CS_SPL_MotionCompensation( ...
  %         reconstructed_frame, reference_frame, block_size,accuracy);
  motion_compensated_image = motion_compensated_image1*weight + motion_compensated_image2*(1-weight);
  %   motion_compensated_image = (motion_compensated_image1+motion_compensated_image2)/2;
  reconstructed_frame = ...
    MC_CS_PL_ResidualReconstruction(motion_compensated_image,...
    y, Phis, num_rows, num_cols);
end