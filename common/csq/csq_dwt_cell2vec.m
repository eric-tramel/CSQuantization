function v = csq_dwt_cell2vec(W)
% A helper function. Converts the cell array DWT to a vector
% The current version uses a slow procedure of expanding v on
% each iteration.

%% Original
% % Get number of levels
% L = length(W);
% v = [];
% for l=1:(L-1)
% 		w = cell2mat(W{l});
% 		v = vertcat(v,w(:));
% end
% v = vertcat(v,W{L}(:));

%% Speed Optimization
% Get number of levels
L = length(W);
sbsize = size(W{1}{1});
imsize = sbsize.*2;
v = zeros(imsize(1)*imsize(2),1);
startidx = 1;
for l=1:(L-1)
        endidx = startidx + 3*sbsize(1)*sbsize(2) - 1;
        w = [W{l}{1} W{l}{2} W{l}{3}];
        v(startidx:endidx) = w(:);
        sbsize = sbsize./2;
        startidx = endidx + 1;
end

v(startidx:end) = W{L}(:);