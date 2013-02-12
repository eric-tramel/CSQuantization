function q = createQuantizer(x, binTable)
% CREATEQUANTIZER Creates quantizer object to be used with state evolution
% or relaxed belief propagation or approximate message passing algorithms.
%
% q = createQuantizer(x, binTable)
%
% Input:
% - x: quantizer decision boundaries.
% - binTable: table of quantizer output binning
%
% e.g. x = [-1e3, -1, 0, 1, 1e3] for practical reasons we replace Inf by 1e3
% e.g. regular N-level quantizer: binTable = repmat((1:N)', 1, 2);
% e.g. L=2 levels per bin with B=2 bins: binTable = [(1:4)', [1;2;1;2]];
%
% Output:
% - q: structure with following fields
% -- q.nLevels: number of levels of underlying regular quantizer 
% (nLevels = size(binTable, 1))
% -- q.nBins: number of bins (nBins = max(binTable(:, 2)))
% -- q.inverse: mapping from index to set of disjoint intervals
% -- q.x: quantizer decision boundaries
%
% Note: q.inverse is a cell array where each component is an array of
% decision boundaries. E.g. inverse{iBin} = [x_0, x_1, x_5, x_6], where
% x_0 < x_1 < x_5 < x_6 are real numbers.
%
% Ulugbek Kamilov, 2010, STIR @ MIT.
%
% See also STATEEVOLUTION, RECONSTRUCTRBP, RECONSTRUCTAMP.

% Number of levels (not-bins)
nLevels = size(binTable, 1);

% Number of bins
nBins = max(binTable(:, 2));

% Initialize the table of inverse mappint  Q^-1(y)
inverse = cell(nBins, 1);

% Go through each level and asign it to a bin
for l = 1:nLevels
    % Extract bin index
    iBin = binTable(l, 2);
    
    % Place in the bin
    inverse{iBin} = [inverse{iBin}, [x(l), x(l+1)]];
end

% Save in the quantizer object
q.nLevels = nLevels;
q.nBins = nBins;
q.inverse = inverse;
q.binTable = binTable;
q.x = x;