function output = get_mean_difference_Bootstrap_t(input_data)
% output = get_mean_difference_Bootstrap_t(input_data)
% custom wrapper function for computing the mean difference between two datasets
% using a provided bootstrap-t implementation.
%   input_data : cell array containing input data matrices
%       input_data{1} - Xinput1
%       input_data{2} - Yinput1
%       input_data{3} - Xinput2
%       input_data{4} - Yinput2
%   output : cell array containing computed results matrices
%       output{1} - mean difference between input_data{1} and input_data{3}
%       output{2} - mean of input_data{1}
%       output{3} - mean of input_data{3}
%
% Giulio Matteucci 2023

% compute the mean of input data 1
mean1 = mean(input_data{1});
% compute the mean of input data 2
mean2 = mean(input_data{3});
% compute the mean difference
if input_data{6}
    mean_diff = mean1 - mean2;
else
    mean_diff = abs(mean1 - mean2);
end
% store results
output{1} = mean_diff;
output{2} = mean1;
output{3} = mean2;

end
