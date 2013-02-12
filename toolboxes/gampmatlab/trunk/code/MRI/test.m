addpath('..\main');
addpath('demo_cs_mri');

% Run the example code up to the point where the parameters are generated
example1_2d_brain_mri_gamp;

% Construct the mriTrans

mriTrans = MriLinTrans(param);
[nout,nin] = mriTrans.size();

% Add this method to extract the components sizes of the outputs
[noutF, noutW, noutTV] = mriTrans.getCompSize();

% Generate a Gaussian input with unit variance
xvar = 1;
x = sqrt(xvar)*randn(nin,1);

% Compute output and predicted variance
y    = mriTrans.mult(x);
yvar = mriTrans.multSq( repmat(xvar,nin,1) );

% Compare measured and predicted variances
% These should be roughly the same.  You can repeat for the wavelet and TV
yvarFtrue = mean( abs(y(1:noutF)).^2 );
yvarF     = mean( yvar(1:noutF) );