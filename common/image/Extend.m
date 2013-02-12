function XE = Extend(X)

[num_rows num_cols] = size(X);
XE = [X(1,:); X; X(num_rows,:)];
XE = [XE(:,1) XE XE(:,num_cols)];
