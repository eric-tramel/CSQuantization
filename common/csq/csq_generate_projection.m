function [A AT] = csq_generate_projection(proj_name,params)
% [A AT] = csq_generate_projection('projection name',params)
%
% Generate forward and transpose function handles for the given
% projection type using the specified parameters. Params must 
% contain all parameters required for the specific projection being
% used as well as the following parameters specific to 
% csq_generate_projection.
%
% Parameter fields:
%
%	block_based:	TRUE, Use block-based CS Acquisition
%			FALSE, Full signal sampling
%	block_dim:	Vector of block size dimensionality.
%	subrate:	Value in [0,1], specify the number of measurements
%			in relation to ambient dimensionality.
%	rand_seed:	Value to set the random number generate (optional).
%	imsize:		Vector of image dimensionality
%	N: 		Dimensionality of signal

%% Variables
block_mode = 0;
rand_seed = java.lang.System.currentTimeMillis; % Current UTC in milliseconds

%% Input checking
csq_required_parameters(params,'subrate');
subrate = params.subrate;

FLAG_Nspecified = 0;
FLAG_imsizespecified = 0;
FLAG_block_dimspecified = 0;

if isfield(params,'N')
	N = params.N;
	FLAG_Nspecified = 1;
end

if isfield(params,'imsize')
	imsize = params.imsize;
	FLAG_imsizespecified = 1;
end

if isfield(params,'block_based')
	block_mode = params.block_based;
end

if isfield(params,'block_dim')
	block_dim = params.block_dim;
	FLAG_block_dimspecified = 1;
end

if isfield(params,'rand_seed')
    rand_seed = params.rand_seed;
end

%% Input handling
% Error producing inputs
if ~FLAG_Nspecified && ~FLAG_imsizespecified
	return_str = sprintf('No signal dimensionality specified.');
	error('csq_generate_projection:NoDim',return_str);
end

if block_mode && ~FLAG_imsizespecified
	return_str = sprintf('Need image dimensions for block acquisition mode.');
	error('csq_generate_projection:NoDim',return_str);
end

if block_mode && ~FLAG_block_dimspecified
	return_str = sprintf('Need block dimensions for block acquisition mode.');
	error('csq_generate_projection:NoDim',return_str);
end

if block_mode && (mod(imsize(1),block_dim(1)) ~= 0 || mod(imsize(2),block_dim(2)) ~= 0) 
	return_str = sprintf('Block dimensions do not divide into image dimensionality evenly.');
	error('csq_generate_projection:DimMismatch',return_str);
end

% Inferring from inputs
if FLAG_imsizespecified && ~FLAG_Nspecified && ~block_mode
	N = imsize(1)*imsize(2);
end

if block_mode
	N = block_dim(1)*block_dim(2);
    Nb = imsize(1)/block_dim(1) * imsize(2)/block_dim(2);
end

%% Handle types of projection
% Get the number of measurements
M = round(subrate*N);

% Set the RNG seed
s = RandStream('mcg16807','Seed',mod(rand_seed,2^32));
RandStream.setDefaultStream(s);

switch proj_name
case 'srm-blk'
	csq_required_parameters(params,'blksize','trans_mode');
	rand_vect = randperm(N)';
	select_vect = randperm(N);
	select_vect = select_vect(1:M);
	Phi   = @(z) blk_f1d(z,select_vect,rand_vect,params.trans_mode,params.blksize);
	PhiT= @(z) blk_t1d(z,N,select_vect,rand_vect,params.trans_mode,params.blksize);

	if block_mode
        A = @(z) vectorize( batch_projection(Phi,...
                                             im2col(reshape(z,imsize),block_dim,'distinct'),...
                                             M,Nb)) ;
        AT = @(z) vectorize(col2im( batch_projection(PhiT,reshape(z,[M Nb]),N,Nb),block_dim,imsize,'distinct'));
	else
		A = Phi;
		AT = PhiT;
	end

case 'srm-fft'
	rand_vect = randperm(N)';
	select_vect = randperm(round(N/2)-1)+1;
	select_vect = select_vect(1:round(M/2))';
	Phi   = @(z) fft1d_f(z, select_vect, rand_vect);
	PhiT = @(z) fft1d_t(z, N, select_vect, rand_vect);

	if block_mode
        A = @(z) vectorize( batch_projection(Phi,...
                                             im2col(reshape(z,imsize),block_dim,'distinct'),...
                                             M,Nb)) ;
        AT = @(z) vectorize(col2im( batch_projection(PhiT,reshape(z,[M Nb]),N,Nb),block_dim,imsize,'distinct'));
	else
		A = Phi;
		AT = PhiT;
	end

case 'gaussian'
    if M <= N
        Phi = orth(randn(N,M))';
        PhiT = Phi';
    else
        Phi = orth(randn(M,N));
        PhiT = Phi';
    end
    
    if block_mode
        A = @(z) vectorize(Phi*im2col(reshape(z,imsize),block_dim,'distinct'));
        AT = @(z) vectorize(col2im(PhiT*reshape(z,[M Nb]),block_dim,imsize,'distinct'));
    else
        A = @(z) Phi*z;
        AT = @(z) PhiT*z;
    end
    
% case 'binary'

otherwise
	return_str = sprintf('Projection "%s" is unsupported.',proj_name);
	error('csq_generate_projection:UnsupportedTransform',return_str);
end	



%----------------------------------------------------
function y = batch_projection(A,x,M,B)
	y = zeros(M,B);
	for i=1:B
		y(:,i) = A(x(:,i));
    end
    
function v = vectorize(y)
	v = y(:);