function v = csq_dwt_cell2vec(W)
% A helper function. Converts the cell array DWT to a vector
% The current version uses a slow procedure of expanding v on
% each iteration.

% Get number of levels
L = length(W);
v = [];
for l=1:(L-1)
		w = cell2mat(W{l});
		v = vertcat(v,w(:));
end

v = vertcat(v,W{L}(:));

