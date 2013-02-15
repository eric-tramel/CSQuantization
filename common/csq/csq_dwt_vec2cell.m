function W = csq_dwt_vec2cell(v,num_rows,num_cols,L)
% A helper function. Converts the given rasterized DWT vector
% to the cell format the WaveletSoftware package wants.

sbsize = [num_rows num_cols];

for l=1:L
	sbsize = sbsize ./ 2;
	sbN = sbsize(1)*sbsize(2);

	W{l}{1} = reshape(v(1:sbN),sbsize);
	v = v(sbN+1:end); % Gobble
	W{l}{2} = reshape(v(1:sbN),sbsize);
	v = v(sbN+1:end); % Gobble
	W{l}{3} = reshape(v(1:sbN),sbsize);
	v = v(sbN+1:end); % Gobble
end

W{l+1} = reshape(v,sbsize);
