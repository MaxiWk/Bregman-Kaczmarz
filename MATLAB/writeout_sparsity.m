% writes sparsity data into a txt file
% requires size(array_on_axis,2) == 1 and 
% size(array_on_axis,1) == size(data,1) 
% header should be a cell array with length(header(2:end)) == size(data, 2)
% and header(1) specifies array_on_xaxis
function writeout_sparsity(writepath, header, min_nnz, median_nnz, max_nnz)   
    writecell(header,  writepath, 'Delimiter', 'tab', 'FileType', 'text');
    writematrix([min_nnz; median_nnz; max_nnz], writepath,...
            'WriteMode', 'append', 'Delimiter', 'tab', 'FileType', 'text');            
end

