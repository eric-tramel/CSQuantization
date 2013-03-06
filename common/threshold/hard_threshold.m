function x_check = hard_threshold(x,lambda,num_rows,num_cols)
x_check = x;
threshold = lambda * sqrt(2 * log(num_rows * num_cols)) * ...
    (median(abs(x(:))) / 0.6745);
x_check(abs(x_check) < threshold) = 0;