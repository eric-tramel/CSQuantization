function W = bivariate_shrinkage(W,lambda,windowsize,end_level)
% thresh = bivariate_shrinkage(W,lambda,windowsize,end_level)
% Perform bivariate shrinkage on the cell-formated wavelet coefficients, W.

% I don't know abou the time issues, but we need to make some 
% allowances for conflicts in Octave with the 'expand' command.
if csq_in_octave()
  local_expand = @(x) wavelet_expand(x);
else
  local_expand = @(x) expand(x);
end

if size(W{1},2) == 3 % for DWT
  windowfilt = ones(1, windowsize)/windowsize;
  
  tmp = W{1}{3};
  Nsig = median(abs(tmp(:)))/0.6745;
  
  for scale = 1:end_level
    for dir = 1:3
      Y_coefficient = W{scale}{dir};
      
      Y_parent = W{scale+1}{dir};
      
      Y_parent = local_expand(Y_parent);
      
      Wsig = conv2(windowfilt, windowfilt, (Y_coefficient).^2, 'same');
      Ssig = sqrt(max(Wsig-Nsig.^2, eps));
      
      T = sqrt(3)*Nsig^2./Ssig;
      
      W{scale}{dir} = bishrink(Y_coefficient, Y_parent, T*lambda);
    end
  end
else %size(W{1},2) == 2 for DDWT
  windowsize  = 3;
  windowfilt = ones(1, windowsize)/windowsize;
  
  tmp = W{1}{1}{1}{1};
  Nsig = median(abs(tmp(:)))/0.6745;
  
  for scale = 1:end_level
    for dir = 1:2
      for dir1 = 1:3
        Y_coef_real = W{scale}{1}{dir}{dir1};
        Y_coef_imag = W{scale}{2}{dir}{dir1};
        Y_parent_real = W{scale+1}{1}{dir}{dir1};
        Y_parent_imag = W{scale+1}{2}{dir}{dir1};
        Y_parent_real  = local_expand(Y_parent_real);
        Y_parent_imag  = local_expand(Y_parent_imag);
        
        Wsig = conv2(windowfilt, windowfilt, (Y_coef_real).^2, 'same');
        Ssig = sqrt(max(Wsig-Nsig.^2, eps));
        
        T = sqrt(3)*Nsig^2./Ssig;
        
        Y_coef = Y_coef_real + sqrt(-1)*Y_coef_imag;
        Y_parent = Y_parent_real + sqrt(-1)*Y_parent_imag;
        Y_coef = bishrink(Y_coef, Y_parent, T*lambda);
        
        W{scale}{1}{dir}{dir1} = real(Y_coef);
        W{scale}{2}{dir}{dir1} = imag(Y_coef);
      end
    end
  end
end