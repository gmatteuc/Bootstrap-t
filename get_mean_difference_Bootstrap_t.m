function output = get_mean_difference_Bootstrap_t(input_data)

%GET_MEAN_DIFFERENCE_BOOTSTRAP_T Custom wrapper function for computing the
% mean difference between two datasets using the provided bootstrap-t
% implementation (get_Bootstrap_t_ci_serial).
%
% INPUTS:
% input_data: cell array, containing the input data matrices
% input_data{1} - Xinput1
% input_data{2} - Yinput1
% input_data{3} - Xinput2
% input_data{4} - Yinput2
% Note: other elements of input_data are not used in this function, they are here 
% just to illustrate the handling of extra inputs by "get_Bootstrap_t_ci" functions
%
% OUTPUTS:
% output: cell array, containing the computed results matrices
% output{1} - mean difference between input_data{1} and input_data{3}
% output{2} - mean of input_data{1}
% output{3} - mean of input_data{3}
%
% USAGE:
% To use this function with the provided bootstrap-t implementation
% (get_Bootstrap_t_ci_serial or get_Bootstrap_t_ci_parallel), pass an handle 
% to it to the bustrap function, e.g. Bfunc = @get_mean_difference_Bootstrap_t
% In order to use the code with your own custom function you can replace the outputs
% with any other calculation over the inputs (then just make sure to
% perform the bootstrap resemplings on the correct dimension of each input
% matrix!)
%
% AUTHOR: [Your name]
% DATE: [Current date]

    % compute the mean of input data 1 - output 2
    mean1 = mean(input_data{1});

    % compute the mean of input data 2 - output 3
    mean2 = mean(input_data{3});
    
    % compute the mean difference - output 1
    if input_data{6}
        mean_diff = mean1 - mean2;
    else
        mean_diff = abs( mean1 - mean2 );
    end
    
    % store results (output 1, 2, 3)
    output{1}=mean_diff;
    output{2}=mean1;
    output{3}=mean2;
    
end