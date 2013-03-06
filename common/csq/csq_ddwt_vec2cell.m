function W = csq_ddwt_vec2cell(v,num_rows,num_cols,L)
% A helper function. Converts the given rasterized DWT vector
% to the cell format the WaveletSoftware package wants.
%   w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..J (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; d2 = 1,2,3 (orientations)
%   w{J+1}{m}{n} - lowpass coefficients
%       d1 = 1,2; d2 = 1,2 

sbsize = [num_rows num_cols];
W = cell(1,L+1);
v = reshape(v, [2*prod(sbsize) 2]);
%% Slow Version
% for l=1:L
% 	sbsize = sbsize ./ 2;
% 	sbN = sbsize(1)*sbsize(2);
% 
% 	W{l}{1} = reshape(v(1:sbN),sbsize);
% 	v = v(sbN+1:end); % Gobble
% 	W{l}{2} = reshape(v(1:sbN),sbsize);
% 	v = v(sbN+1:end); % Gobble
% 	W{l}{3} = reshape(v(1:sbN),sbsize);
% 	v = v(sbN+1:end); % Gobble
% end
% 
% W{l+1} = reshape(v,sbsize);

%% Fast Version?
for l=1:L

  sbsize = sbsize ./ 2;
  sbN = prod(sbsize);
  i = 1;
  W{l}{i}{1}{1} = reshape(v(1:sbN,i),sbsize);
  W{l}{i}{1}{2} = reshape(v(sbN+1:2*sbN,i),sbsize);
  W{l}{i}{1}{3} = reshape(v(2*sbN+1:3*sbN,i),sbsize);
  
  W{l}{i}{2}{1} = reshape(v(3*sbN+1:4*sbN,i),sbsize);
  W{l}{i}{2}{2} = reshape(v(4*sbN+1:5*sbN,i),sbsize);
  W{l}{i}{2}{3} = reshape(v(5*sbN+1:6*sbN,i),sbsize);  
  
  i = 2;
  W{l}{i}{1}{1} = reshape(v(1:sbN,i),sbsize);
  W{l}{i}{1}{2} = reshape(v(sbN+1:2*sbN,i),sbsize);
  W{l}{i}{1}{3} = reshape(v(2*sbN+1:3*sbN,i),sbsize);
  
  W{l}{i}{2}{1} = reshape(v(3*sbN+1:4*sbN,i),sbsize);
  W{l}{i}{2}{2} = reshape(v(4*sbN+1:5*sbN,i),sbsize);
  W{l}{i}{2}{3} = reshape(v(5*sbN+1:6*sbN,i),sbsize);  
 
  v = v(6*sbN+1:end,:); % Gobble
end

W{l+1}{1}{1} = reshape(v(1:sbN,1),sbsize);
W{l+1}{1}{2} = reshape(v(sbN+1:2*sbN,1),sbsize);
W{l+1}{2}{1} = reshape(v(1:sbN,2),sbsize);
W{l+1}{2}{2} = reshape(v(sbN+1:2*sbN,2),sbsize);
