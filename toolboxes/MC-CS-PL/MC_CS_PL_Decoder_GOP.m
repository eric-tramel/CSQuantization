function [reconstructed_frames reconstructed_frames_indep] = ...
  MC_CS_PL_Decoder_GOP(y, Phis, GOP_size, num_rows, num_cols)
% global frames
Phi_key = Phis{1};
Phi_key_t = Phis{2};
Phi = Phis{3};
Phi_t = Phis{4};

ind_center_frame = GOP_size/2 + 1;

disp('< 2D-CS-PL Initial Reconstruction >');

recon_2d = cell(1,GOP_size+1);
recon_2d{1} = CS_PL_DDWT_Decoder(y{1}, Phi_key, Phi_key_t, num_rows, ...
  num_cols);
for i = 2:GOP_size
  recon_2d{i} = CS_PL_DDWT_Decoder(y{i}, Phi, Phi_t, num_rows, ...
    num_cols);
end
recon_2d{GOP_size+1}= CS_PL_DDWT_Decoder(y{end}, Phi_key, Phi_key_t, ...
  num_rows, num_cols);

% for i =1:GOP_size + 1
%   psnrs_2d(i) = PSNR(recon_2d{i},frames{i});
% end
% psnrs_2d
% mean(psnrs_2d)

reconstructed_frames_indep = recon_2d(1:end-1);

disp('< Forward Processing >');
[recon_forward_frames recon_2d] = ...
  MC_CS_PL_Decoder_ForwardBackward(y(1:ind_center_frame-1), Phis, ...
  num_rows, num_cols, recon_2d(2:ind_center_frame-1), recon_2d(1:end));

disp('< Backward Processing >' );
[recon_backward_frames recon_2d]= ...
  MC_CS_PL_Decoder_ForwardBackward( y(GOP_size+1:-1:ind_center_frame+1), ...
  Phis, num_rows, num_cols, recon_2d(GOP_size:-1:ind_center_frame+1), recon_2d(end:-1:1));

disp('< Bidirectional Processing for Center frame >');
recon_center_frame = ...
  MC_CS_PL_Decoder_Center(y{ind_center_frame}, ...
  recon_forward_frames{end}, recon_backward_frames{end}, Phis, recon_2d{ind_center_frame}, recon_2d);

reconstructed_frames = ...
  [recon_forward_frames recon_center_frame recon_backward_frames(end:-1:2)];
