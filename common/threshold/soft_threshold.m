function x_check = soft_threshold(x,lambda,num_rows,num_cols)
threshold = lambda * sqrt(2 * log(num_rows * num_cols)) * ...
    (median(abs(x(:))) / 0.6745);
tmp=abs(x)-threshold;
x_check = sign(x).*tmp.*(tmp>0);
  
