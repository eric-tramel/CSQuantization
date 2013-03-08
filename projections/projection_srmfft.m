function [A AT B BT] = projection_srmfft(M,N,is_blocked,imsize,block_dim)
% BROKEN -- Until further notice, use only projection_srmfft
A  = [];
AT = [];
B  = [];
BT = [];
disp('projection_srmfft: Broken, do not use me');


%% Archived Code
% 	if nargin < 3
% 		is_blocked = 0;
% 		imsize = [];
% 		block_dim = [];
% 		Nb = 0;
% 	end
% 	
% 	% Split up projection sets to allow for Subrates > 1
% 	Mleft = M;
% 	i = 1;
% 	while Mleft >= N
% 		thisM = (N - 1);
% 		[B{i} BT{i}] = single_projector_set(thisM,N,is_blocked,imsize,block_dim,Nb);
% 		i = i + 1;
% 		Mleft = Mleft - thisM;
% 	end
% 
% 	% If there are left over measurements
% 	if Mleft ~= 0
% 		[B{i} BT{i}] = single_projector_set(Mleft,N,is_blocked,imsize,block_dim,Nb);
% 	end
% 
% 	% Set up forward transform
% 	A = @(z) srm_forward_batch(B,z);
% 
% 	% Set up transpose transform, need to rearrange y's to account
% 	% To do this, we need to sum together all the individual projections
%     if is_blocked
%         AT = @(z) srm_transpose_batch(BT,z,N*Nb,Nb,i,M./N);
%     else
%         AT = @(z) srm_transpose_batch(BT,z,N,1,i,M./N);
%     end
% 
% 
% %----------------------------------------------------
% function [A AT] = single_projector_set(M,N,is_blocked,imsize,block_dim,Nb)
%         rand_vect = randperm(N)';
%         select_vect = randperm(round(N/2)-1)+1;
%         select_vect = select_vect(1:floor(M/2))';
%         Phi   = @(z) fft1d_f(z, select_vect, rand_vect);
%         PhiT = @(z) fft1d_t(z, N, select_vect, rand_vect);
% 
% 	if is_blocked
%         A = @(z) srm_forward(Phi,z,M,imsize,block_dim,Nb);
%         AT = @(z) srm_transpose(PhiT,z,M,N,block_dim,imsize,Nb);
%     else
% 		A = Phi;
% 		AT = PhiT;
%     end
% 
% 
% function y = srm_forward(Phi,x,M,imsize,block_dim,Nb)
%     y = reshape(x,imsize);
%     y = im2col(y,block_dim,'distinct');
%     y = projection_batch(Phi,y,M,Nb);
%     y = y(:);
% 
% function x = srm_transpose(PhiT,y,M,N,block_dim,imsize,Nb)
%     x = reshape(y,[M,Nb]);
%     x = projection_batch(PhiT,x,N,Nb);
%     x = col2im(x,block_dim,imsize,'distinct');
%     x = x(:);
%     
% % function y = batch_projection(A,x,M,B)
% % % 	y = zeros(M,B);
% %     for i=1:B
% % 		y(:,i) = A(x(:,i));
% %     end
% 
% function x = srm_transpose_batch(AT,y,N,Nb,num_projs,rescale)
% 	x = zeros(N,1);
% 	i = 1;
% 	while i < num_projs
% 		x = x + AT{i}(y(1:(N-Nb)));
% 		y = y((N-Nb+1):end);
% 		i = i + 1;
% 	end
% 	x = x + AT{i}(y);
%     x = x ./ rescale;
% 
% function y = srm_forward_batch(A,x)
%     % Assuming A is a cell array
%     y = [];
%     for i=1:length(A)
%         y = vertcat(y,A{i}(x));
%     end