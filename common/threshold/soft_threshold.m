function thresh = soft_threshold(x,t)
  tmp=abs(x)-t;
  thresh = sign(x).*tmp.*(tmp>0);