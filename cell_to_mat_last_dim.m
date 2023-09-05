function output_mat = cell_to_mat_last_dim(input_cell)
% cell_to_mat_last_dim(input_cell)
% convert a cell array to a matrix by concatenating along an extra dimension.
%   input_cell : cell array containing matrices to concatenate
%   output_mat : concatenated matrix along the last dimension
%
% Giulio Matteucci 2023

% initialize output matrix
output_mat = [];
% get the number of rows in the cell array
num_rows = size(input_cell, 1);
% loop through each cell to concatenate matrices along the last dimension
for i = 1:num_rows
    current_matrix = input_cell{i};
    num_dims = ndims(current_matrix);
    output_mat = cat(num_dims + 1, output_mat, current_matrix);
end

end
