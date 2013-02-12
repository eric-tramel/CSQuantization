function X = loadImageBatch(filenames,color_fmt)
% X = loadImageBatch(filenames)
% Takes a cell array of image file names and loads all images into the cell
% array, X.

num_files = length(filenames);
X = cell(num_files,1);

for i=1:num_files
    X{i} = double(imread(filenames{i}));
    color_dim = size(X{i},3);
    
    if color_dim == 3
    % Operations on three color dimensions
        switch color_fmt
            case 'lum'
                X{i} = rgb2yuv(X{i});
                X{i} = squeeze(X{i}(:,:,1));
            case 'yuv'
                X{i} = rgb2yuv(X{i});
        end
    else
    % Operations when only one color dimension is read
    % from the image.
        switch color_fmt
            case 'rgb'
                X{i} = lum2rgb(X{i});
        end
    end
    
    
end