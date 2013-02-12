%
% function [reconstructed_current_frame recon_indep] = ...
%     MC_BCS_SPL_Decoder(reference_frame, y, Phi, num_rows, num_cols)
%
%   This function performs MC-BCS-SPL method to reconstruct current frame
%   based on reference frame. Multihypothesis initialization of independent
%   reconstruction and residual reconstruction without motino is utilized
%   for the initial reconstruction of current frame to obtain motion
%   vectors. After that, residual reconstruction is performed which enables
%   better reconstruction than direct reconstruction because it is sparser.
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

function reconstructed_current_frame = ...
  MC_CS_PL_Decoder(reference_frame, recon_3d, y, Phis, num_rows, num_cols, recon_3ds)

subrate = Phis{6};

weight = subrate*2 - 0.2;
% weight = 0.5;
weight(weight>0.8) = 0.8;
weight(weight<0.2) = 0.2;

recon_without_motion = ...
  MC_CS_PL_ResidualReconstruction(reference_frame, y, Phis, ...
  num_rows, num_cols);
reconstructed_current_frame = recon_3d*weight + recon_without_motion*(1-weight);

accuracy = 'fp';
num_iter = 1;
block_size = 32;
window_size = 7;

recon_3ds = [recon_3ds {reconstructed_current_frame}];
for i = 1:num_iter
  
  motion_compensated_frame1=  memc_telescope(recon_3ds, reconstructed_current_frame, block_size, window_size);
  %   motion_compensated_frame1=  motion_estimation_compensation(recon_3ds{2}, reconstructed_current_frame, block_size, window_size);
  %   motion_compensated_frame2 = optical_flow(reconstructed_current_frame, reference_frame);
  motion_compensated_frame2=  motion_estimation_compensation(reference_frame, reconstructed_current_frame, block_size, window_size);
  %   motion_compensated_frame2 = MC_CS_SPL_MotionCompensation( ...
  %     reconstructed_current_frame, reference_frame, block_size,accuracy);
  motion_compensated_frame = motion_compensated_frame1*weight + motion_compensated_frame2*(1-weight);
  reconstructed_current_frame = MC_CS_PL_ResidualReconstruction( ...
    motion_compensated_frame, y, Phis, num_rows, num_cols);
  
  
end
