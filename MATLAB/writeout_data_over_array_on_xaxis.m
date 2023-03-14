% writes data into a txt file
% requires size(array_on_axis,2) == 1 and 
% size(array_on_axis,1) == size(data,1) 
% header should be a cell array with length(header(2:end)) == size(data, 2)
% and header(1) specifies array_on_xaxis
function writeout_data_over_array_on_xaxis(writepath, header, array_on_xaxis, data) 

    if size(array_on_xaxis, 1) == 1
       array_on_xaxis = reshape(array_on_xaxis, size(array_on_xaxis, 2), size(array_on_xaxis, 1));
    end
    
    xy = [array_on_xaxis, data];
            writecell(header,  writepath, 'Delimiter', 'tab', 'FileType', 'text');
            writematrix(xy, writepath,...
                    'WriteMode', 'append', 'Delimiter', 'tab', 'FileType', 'text');     
end

