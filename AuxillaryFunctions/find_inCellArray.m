function [cell_Adresses, cellIn_Adresses] = find_inCellArray(cellArray, opt_condition)
% [cellAdress, cellInsideAdress] = find_inCellArray(cellArray, opt_condition)
% Similar to MATLAB built-in function 'find' but for cell arrays.
% The main differnce is the contion should be passed to this function as
% 'function handel'.
% 
% EXAMPLE:
% a = 2;
% [cell_Adresses, cellIn_Adresses] = find_inCellArray(cellArray, @(X)  (x == a))
% ------
% Input:
% 1      cellArray: cell array :D
%        The cell array that you gonna search in it
% 2      opt_condition: function handel
%        EXP: @isnan, @(x) x == 2
%        This function handle express the condition that you'd like to know
%        which component of your cell array is satisfying this condition.
%        This operator will apply in each individual component of the cell
%        arrray and return they logical indices.
% 
% Output:
% 1     cell_Adresses: nCF*nDim
%       Addresses of cell component who satisfied within in the condition
%       (expressed in 'opt_condition'). The addresses will be stored in
%       matrix in which each row provide the addresses of the specific cell
%       component that satisfied the condition. 
%       The number of columns depends on the size of the 'cellArray'. If the
%       cellArray is n-D you will have nDim columns.
%       EXP: if the input cellArray was {3*4} and 5 component of it satisfy the condition then cell_Adresses will be 5*2                      
% 1     cellIn_Adresses: logical matrix 
% ------
% see also FIND
% ------
% potential improvments:
% (1) how you would like format the output (logical incies, indices,
% row/col), format similar to find
% (2) what about if you have nested cell within array
% ------
% Code Info:
%   creation: 2015-06-23 by ShS -> shervin.safavi@gmail.com
%   modification: 

nCellComponenet = numel(cellArray);

% go through all the cell components and check if the condition implemented
% in opt_condition is satisfied 
% and
% store the indices based on what is stored in 'logicalIndices'
nC_sf = 0; % to count ceontent of how many cells satisfy the condition
for iC = 1 : nCellComponenet
    logicalIndex = opt_condition(cellArray{iC}); % check the condition within each cell   
    if sum(logicalIndex(:)) >= 1
        nC_sf = nC_sf + 1; % nCell_satisfiedCondition
        cellIn_logicalIndex{nC_sf,1} = logicalIndex;
        cell_linearIndex(nC_sf) = iC;
%         cell_Adresses(nC_sf, :) = ind2sub(size(cellArray), iC);
        
    end
end

% here you can think about how you gonna send it out (e.g. linear index, logical index)
cellIn_Adresses = cellIn_logicalIndex;
[cell_Adresses{1:ndims(cellArray),1}] = ind2sub(size(cellArray), cell_linearIndex);
cell_Adresses = cell2mat(cell_Adresses)';



% Why cell fun was not used?
% http://stackoverflow.com/questions/18284027/cellfun-versus-simple-matlab-loop-performance
% for custom functions is slower than normal for loops
% that's why not:
% logicalIndices = cellfun(opt_condition, cellArray, 'UniformOutput', false); 