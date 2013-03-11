% projection_example.m
%
% A small demonstration of how to use the projection tools. This script
% also serves to test the projection tools.
clear

% Experiment variables
proj_type = 'gaussian';
block_based = 1;
S1 = 0.25;
S2 = 0.5;
S3 = 2;


% Dependency Loading
csq_deps('srm','proj','common');

% Data load
X = csq_load_data('image','lena.jpg');
x = X(:);

% Common projection parameters
params.block_based = block_based;
params.block_dim = [32 32];
params.trans_mode = 'BWHT';
params.blksize = 32;
params.rand_seed = 1;
params.imsize = size(X);
params.N = length(x);

% Specific projection parameters
params1 = params;
params1.subrate = S1;

params2 = params;
params2.subrate = S2;

params3 = params;
params3.subrate = S3;

% Generate projection sets
[A1 AT1] = csq_generate_projection(proj_type,params1);
[A2 AT2] = csq_generate_projection(proj_type,params2);
[A3 AT3] = csq_generate_projection(proj_type,params3);

% Get transpose images
X1 = reshape(AT1(A1(x)),params.imsize);
X2 = reshape(AT2(A2(x)),params.imsize);
X3 = reshape(AT3(A3(x)),params.imsize);

% Display Results
figure(1);
subplot(1,3,1);
    imagesc(X1); axis image;
    title(sprintf('%s, S = %0.2f',proj_type,S1));
subplot(1,3,2);   
    imagesc(X2); axis image;
    title(sprintf('%s, S = %0.2f',proj_type,S2));
subplot(1,3,3);
    imagesc(X3); axis image;
    title(sprintf('%s, S = %0.2f',proj_type,S3));
colormap(gray);



