
function block = extract_block(image, row, col, block_size)

[num_rows num_cols] = size(image);

block = zeros(block_size, block_size);

if ((row < 1) || (col < 1) || ...
      ((row + block_size - 1) > num_rows) || ...
      ((col + block_size - 1) > num_cols))
  row3 = row;
  for row1 = 1:block_size
    col3 = col;
    for col1 = 1:block_size
      if (row3 < 1)
        row2 = 1;
      else
        if (row3 > num_rows)
          row2 = num_rows;
        else
          row2 = row3;
        end
      end
      if (col3 < 1)
        col2 = 1;
      else
        if (col3 > num_cols)
          col2 = num_cols;
        else
          col2 = col3;
        end
      end
      block(row1, col1) = image(row2, col2);
      col3 = col3 + 1;
    end
    row3 = row3 + 1;
  end
else  
  block = image(row:(row + block_size - 1), col:(col + block_size - ...
      1));
end
          
