function z = RasterScan(X,style)
% This function will apply some kind of raster scanning method on the two
% dimensional image/dataset, X.
[num_rows num_cols] = size(X);

if min([num_rows num_cols]) == 1
    fprintf('[Error] RasterScan(): Input X must be two dimensional.\n');
end

switch style
    case 'columns'
        z = X(:);
    case 'rows'
        X = X';
        z = X(:);
    case 'downup'
        z = down_and_up(X);
    case 'updown'
        z = up_and_down(X);
    case 'rightleft'
        X = X';
        z = down_and_up(X);
    case 'leftright'
        X = X';
        z = up_and_down(X);
    case 'zigzag'
        z = zigzag(X);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = down_and_up(X)
num_cols = size(X,2);
evens = 2:2:num_cols;
X(:,evens) = flipud(X(:,evens));
z = X(:);

function z = up_and_down(X)
num_cols = size(X,2);
odds = 1:2:num_cols;
X(:,odds) = flipud(X(:,odds));
z = X(:);

function z = zigzag(X)
[num_rows num_cols] = size(X);
X = fliplr(X);
z = [];
for i=(num_cols-1):-1:-(num_cols-1)
    line = diag(X,i);
    if mod(i,2)
        line = flipud(line(:));
    end
    z = [z;line];
end