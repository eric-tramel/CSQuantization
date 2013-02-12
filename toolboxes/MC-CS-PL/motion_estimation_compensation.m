function  motion_compensated_image= motion_estimation_compensation(reference_image, current_image, B, W)

[num_rows num_cols] = size(current_image);
block_num_rows = floor(num_rows / B);

reference_image = symextend(reference_image,2*W);
current_image = im2col(current_image,[B B],'distinct');

for block_row = 1:B:num_rows
  for block_col = 1:B:num_cols
  
  index = ((block_col-1)/B)*block_num_rows + (block_row-1)/B + 1;
  current_block = current_image(:,index);  
  
  % i,j center of the window
  i = block_row + 2*W;
  j = block_col + 2*W;
  
  % extract all blocks from windows
  window_ref = reference_image(i-W:i+W+B-1,j-W:j+W+B-1);
  H = im2col(window_ref,[B B],'sliding');
  % obtaining mean absolute error for all blocks extracted from the window
  MAEs = mean(abs(bsxfun(@minus, H, current_block)));
  [val ind] = min(MAEs);
  
%     W_size = 2*W+1;
%   [temp_row temp_col] = ind2sub([W_size W_size],ind);

%   temp_col2 = temp_col + j - W -1;
%   temp_row2 = temp_row + i - W -1;

%   motion_compensated_image(block_row:block_row+B-1,block_col:block_col+B-1) = ...
%     reference_image(temp_row2:temp_row2+B-1,temp_col2:temp_col2+B-1);

motion_compensated_image(block_row:block_row+B-1,block_col:block_col+B-1) = ...
     reshape(H(:,ind),[B B]);
  
  end
end
