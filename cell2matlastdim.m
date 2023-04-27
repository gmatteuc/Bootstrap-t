function outputMat = cell2matlastdim(inputCell)
% convert cell array to matrix concatenating on an extra dimension
N = size(inputCell, 1);
outputMat = [];
for i = 1:N
    currentMatrix = inputCell{i};
    Ndim = ndims(currentMatrix);
    outputMat = cat(Ndim+1, outputMat, currentMatrix);
end
end

