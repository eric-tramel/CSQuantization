function  motion_compensated_image= memc_telescope(reference_frames, current_image, B, W)

[num_rows num_cols] = size(current_image);
block_num_rows = floor(num_rows / B);

ref_size = length(reference_frames);
for i = 1:ref_size
  reference_frames{i} = symextend(reference_frames{i},2*W);
end
current_image = im2col(current_image,[B B],'distinct');

for block_row = 1:B:num_rows
  for block_col = 1:B:num_cols
  
  index = ((block_col-1)/B)*block_num_rows + (block_row-1)/B + 1;
  current_block = current_image(:,index);  
  
  % i,j center of the window
  i = block_row + 2*W;
  j = block_col + 2*W;
  
  % extract all blocks from windows
  H = cell(1,ref_size);
  for k = 1:ref_size
    window_ref = reference_frames{k}(i-W:i+W+B-1,j-W:j+W+B-1);
    H{k} = im2col(window_ref,[B B],'sliding');
  end
  H = cell2mat(H);
  % obtaining mean absolute error for all blocks extracted from the window
  MAEs = mean(abs(bsxfun(@minus, H, current_block)));
  [val ind] = min(MAEs);

  motion_compensated_image(block_row:block_row+B-1,block_col:block_col+B-1) = ...
    reshape(H(:,ind),[B B]);
  
  end
end
