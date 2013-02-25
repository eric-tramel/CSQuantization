function [A AT B BT] = projection_srmblk(M,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb)
	if nargin < 5
		is_blocked = 0;
		imsize = [];
		block_dim = [];
		Nb = 0;
	end
	
	% Split up projection sets to allow for Subrates > 1
	Mleft = M;
	i = 1;
	while Mleft >= N
		thisM = (N - 1);
		[B{i} BT{i}] = single_projector_set(thisM,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb);
		i = i + 1;
		Mleft = Mleft - thisM;
	end

	% If there are left over measurements
	if Mleft ~= 0
		[B{i} BT{i}] = single_projector_set(Mleft,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb);
	end

	% Set up forward transform
	A = @(z) srm_forward_batch(B,z);

	% Set up transpose transform, need to rearrange y's to account
	% To do this, we need to sum together all the individual projections
	AT = @(z) srm_transpose_batch(BT,z,M*Nb,N*Nb,i);


%----------------------------------------------------
function [A AT] = single_projector_set(M,N,trans_mode,blksize,is_blocked,imsize,block_dim,Nb)
	rand_vect = randperm(N)';
	select_vect = randperm(N);
	select_vect = select_vect(1:M);
	Phi   = @(z) blk_f1d(z,select_vect,rand_vect,trans_mode,blksize);
	PhiT  = @(z) blk_t1d(z,N,select_vect,rand_vect,trans_mode,blksize);

	if is_blocked
	    A = @(z) vectorize( batch_projection(Phi,...
	                                 im2col(reshape(z,imsize),block_dim,'distinct'),...
	                                 M,Nb)) ;
	    AT = @(z) vectorize(col2im( batch_projection(PhiT,reshape(z,[M Nb]),N,Nb),block_dim,imsize,'distinct'));
	else
		A = Phi;
		AT = PhiT;
	end


function y = batch_projection(A,x,M,B)
	y = zeros(M,B);
	for i=1:B
		y(:,i) = A(x(:,i));
    end

function x = srm_transpose_batch(AT,y,M,N,num_projs)
	x = zeros(N,1);
	i = 1;
	while i < num_projs
		x = x + AT{i}(y(1:(N-1)));
		y = y(N:end);
		i = i + 1;
	end
	x = x + AT{i}(y);

function y = srm_forward_batch(A,x)
    % Assuming A is a cell array
    y = [];
    for i=1:length(A)
        y = vertcat(y,A{i}(x));
    end

function v = vectorize(y)
	v = y(:);