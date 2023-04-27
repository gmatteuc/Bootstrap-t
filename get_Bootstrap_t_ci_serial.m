function [estimate,estimate_lCI,estimate_uCI,estimate_lSE,estimate_uSE] =...
    get_Bootstrap_t_ci_serial(Bfunc,Binp,Bdim,Brdim,confidence,B,N,seed,printeveryn)

%GET_BOOTSTRAP_T_CI_SERIAL Perform Bootstrap-t (serial) for the estimation of confidence intervals
% and standard error intervals of a summary statistic (any custom function Bfunc).
%
% INPUTS:
% Bfunc: function handle of summary statistic (any custom function Bfunc)
% Binp: cell array, inputs for the Bfunc
% Bdim: vector, dimensions to resample for each input in Binp
% Brdim: vector, resampling dimensions for each input matrix contained in each element of Binp
% confidence: scalar, confidence level (e.g., 0.95 for 95% CI)
% B: integer, number of outer bootstrap repetitions
% N: integer, number of inner bootstrap repetitions
% seed: integer, seed for random number generator
% printeveryn: integer, every how any bootstrap repetitions to print advancement
%
% OUTPUTS:
% estimate: cell array, point estimates of the summary statistic
% estimate_lCI: cell array, lower bound of the confidence interval
% estimate_uCI: cell array, upper bound of the confidence interval
% estimate_lSE: cell array, lower bound of the standard error interval
% estimate_uSE: cell array, upper bound of the standard error interval
%
% EXAMPLE:
% confidence = 0.95;
% B = 2500;
% N = 25;
% Bfunc = @get_mean_difference;
% Binp = {Xinput1, Yinput1, Xinput2, Yinput2, quadrant, takeabs};
% Bdim = [1,1,2,2,0,0]; ---> two independent resampling on X1,Y1 and X2,Y2
% Brdim = [2,2,2,2,0,0]; ---> resampling done on dimension 2 of input matrices X1,Y1,X2,Y2
% [estimate, estimate_lCI, estimate_uCI, estimate_lSE, estimate_uSE] =...
% get_Bootstrap_t_ci_serial(Bfunc,Binp,Bdim,Brdim,confidence,B,N,seed);
%
% AUTHOR: Giulio Matteucci
% DATE: 21/04/2023

%% perform Bootstrap-t (serial)

% get number of outputs
O=length(Bfunc(Binp));

% initialize bootstrap Z variables
BZvals=cell(O,B);
Binp_bo=cell(B,numel(Binp));
Binp_bi=cell(N,numel(Binp));

% initialize estimate structures
estimate=cell(O,1);
estimate_lCI=cell(O,1);
estimate_uCI=cell(O,1);
estimate_lSE=cell(O,1);
estimate_uSE=cell(O,1);

% get observed summary statistics value
Bout=Bfunc(Binp(:));

% nested bootstrap - outer loop ---------------------------------------
Bout_bo=cell(O,B);
for bo_idx=1:B
    tic
    
    % set the seed for the current iteration
    rng(seed + bo_idx);
    
    % get list of resamplings to be done
    requested_resamplings=unique(Bdim);
    % reject "null" resampling label
    requested_resamplings=requested_resamplings(not(requested_resamplings==0));
    % initialize current outer bootstrap run input as default input
    Binp_bo(bo_idx,:)=Binp;
    
    % loop over resamplings
    for requested_resampling_idx=1:numel(requested_resamplings)
        
        % get current resampling label
        current_requested_resampling=requested_resamplings(requested_resampling_idx);
        % get corresponding input structure indeces
        inputs_to_resample=find(Bdim==current_requested_resampling);
        % get current dim to resample
        current_dim_to_resample=Brdim(inputs_to_resample(1));
        % get elements to sample
        current_sample_idx=randsample(size(Binp{inputs_to_resample(1)},current_dim_to_resample),size(Binp{inputs_to_resample(1)},current_dim_to_resample),1);
        
        % loop input structure indeces (i.e. variables to resample together)
        for input_to_resample_idx=1:numel(inputs_to_resample)
            % get current input structure index
            current_input_to_resample=inputs_to_resample(input_to_resample_idx);
            % initialize an index cell array
            current_indexCell = repmat({':'}, 1, ndims(Binp{current_input_to_resample}));
            % replace the specified resampling dimension indices with the elements to sample
            current_indexCell{current_dim_to_resample} = current_sample_idx;
            % get current sample
            current_sample=Binp{current_input_to_resample}(current_indexCell{:});
            % store current sample
            Binp_bo{bo_idx,current_input_to_resample}=current_sample;
        end
        
    end
    
    % compute summary statistics of interest
    Bout_bo(:,bo_idx)=Bfunc(Binp_bo(bo_idx,:));
    
    % nested bootstrap - inner loop  ---------------------------------------
    Bout_bi=cell(O,N);
    for bi_idx=1:N
        
        % set the seed for the current iteration
        rng(seed + bo_idx + bi_idx);
        
        % initialize current inner bootstrap run input as current outer bootstrap run input
        current_Binp_bo=Binp_bo(bo_idx,:);
        Binp_bi(bi_idx,:)=current_Binp_bo;
        
        % loop over resamplings
        for requested_resampling_idx=1:numel(requested_resamplings)
            
            % get current resampling label
            current_requested_resampling=requested_resamplings(requested_resampling_idx);
            % get corresponding input structure indeces
            inputs_to_resample=find(Bdim==current_requested_resampling);
            % get current dim to resample
            current_dim_to_resample=Brdim(inputs_to_resample(1));
            % get elements to sample
            current_sample_idx=randsample(size(current_Binp_bo{inputs_to_resample(1)},current_dim_to_resample),size(current_Binp_bo{inputs_to_resample(1)},current_dim_to_resample),1);
            
            % loop input structure indeces (i.e. variables to resample together)
            for input_to_resample_idx=1:numel(inputs_to_resample)
                % get current input structure index
                current_input_to_resample=inputs_to_resample(input_to_resample_idx);
                % initialize an index cell array
                current_indexCell = repmat({':'}, 1, ndims(Binp{current_input_to_resample}));
                % replace the specified resampling dimension indices with the elements to sample
                current_indexCell{current_dim_to_resample} = current_sample_idx;
                % get current sample
                current_sample=current_Binp_bo{current_input_to_resample}(current_indexCell{:});
                % store current sample
                Binp_bi{bi_idx,current_input_to_resample}=current_sample;
            end
            
        end
        
        % compute summary statistics of interest
        Bout_bi(:,bi_idx)=Bfunc(Binp_bi(bi_idx,:));
        
    end
    % end of inner loop -----------------------------------------------
    
    % compute bootstrap Z
    for output_idx=1:O
        Znum=Bout{output_idx}-Bout_bo{output_idx,bo_idx};
        catmat = cell2matlastdim(Bout_bi(output_idx,:)');
        Zdenom=nanstd(catmat,[],ndims(catmat));
        BZvals{output_idx,bo_idx}=(Znum)./(Zdenom);
    end
    
    % print advancement every printeveryn iterations
    if mod(bo_idx,printeveryn)==0
        toc
        fprintf(['\nbootstrap iteration n=',num2str(bo_idx),' completed \n'])
    end
    
end

for output_idx=1:O
    % get confidence interval
    catmat1 = cell2matlastdim(BZvals(output_idx,:)');
    catmat2 = cell2matlastdim(Bout_bo(output_idx,:)');
    estimate_lCI{output_idx}=Bout{output_idx}-quantile(catmat1,1-(1-confidence)/2,ndims(catmat1)).*std(catmat2,[],ndims(catmat2));
    estimate_uCI{output_idx}=Bout{output_idx}-quantile(catmat1,(1-confidence)/2,ndims(catmat1)).*std(catmat2,[],ndims(catmat2));
    % get standard error interval
    estimate_lSE{output_idx}=Bout{output_idx}-quantile(cell2mat(BZvals(output_idx,:)'),1-(1-0.682)/2).*std(cell2mat(Bout_bo(output_idx,:)'));
    estimate_uSE{output_idx}=Bout{output_idx}-quantile(cell2mat(BZvals(output_idx,:)'),(1-0.682)/2).*std(cell2mat(Bout_bo(output_idx,:)'));
    % get estimate
    estimate{output_idx}=Bout{output_idx};
end

end