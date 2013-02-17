function W = bivariate_shrinkage(W,lambda,windowsize,end_level)
% thresh = bivariate_shrinkage(W,lambda,windowsize,end_level)
% Perform bivariate shrinkage on the cell-formated wavelet coefficients, W.

windowfilt = ones(1, windowsize)/windowsize;

tmp = W{1}{3};
Nsig = median(abs(tmp(:)))/0.6745;

for scale = 1:end_level
  for dir = 1:3
    Y_coefficient = W{scale}{dir};
    
    Y_parent = W{scale+1}{dir};
    
    Y_parent = expand(Y_parent);
    
    Wsig = conv2(windowfilt, windowfilt, (Y_coefficient).^2, 'same');
    Ssig = sqrt(max(Wsig-Nsig.^2, eps));
    
    T = sqrt(3)*Nsig^2./Ssig;
    
    W{scale}{dir} = bishrink(Y_coefficient, Y_parent, T*lambda);
  end
end