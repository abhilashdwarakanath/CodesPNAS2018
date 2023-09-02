parfor i = 1:params.datasets %ndatasets
    [coheron{i}, coheroff{i}] = lfpcoherence(params,lfpall{i},i);
end