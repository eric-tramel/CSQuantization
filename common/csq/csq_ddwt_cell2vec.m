function v = csq_ddwt_cell2vec(W)
% A helper function. Converts the cell array DWT to a vector
% The current version uses a slow procedure of expanding v on
% each iteration.
%   w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..J (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; d2 = 1,2,3 (orientations)
%   w{J+1}{m}{n} - lowpass coefficients
%       d1 = 1,2; d2 = 1,2 

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
sbsize = size(W{1}{1}{1}{1});
imsize = sbsize.*2;
v = zeros(2*prod(imsize),2);
startidx = 1;
for l=1:(L-1)
        endidx = startidx + 6*prod(sbsize) - 1;
        i = 1;
        w = [W{l}{i}{1}{1} W{l}{i}{1}{2} W{l}{i}{1}{3} ...
             W{l}{i}{2}{1} W{l}{i}{2}{2} W{l}{i}{2}{3}];
        v(startidx:endidx,i) = w(:);
        i = 2;
        w = [W{l}{i}{1}{1} W{l}{i}{1}{2} W{l}{i}{1}{3} ...
             W{l}{i}{2}{1} W{l}{i}{2}{2} W{l}{i}{2}{3}];
        v(startidx:endidx,i) = w(:);        
        sbsize = sbsize./2;
        startidx = endidx + 1;
end
i = 1;
w = [W{L}{i}{1} W{L}{i}{2}];
v(startidx:end,i) = w(:);
i = 2;
w = [W{L}{i}{1} W{L}{i}{2}];
v(startidx:end,i) = w(:);

v = v(:);