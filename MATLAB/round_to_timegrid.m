% input: data matrix of dimension "(floor(maxiter)/iter_save) x ..."
%        runtime_over_iter: 1d array with saved runtime for each iteration,
%                           should be monotonely increasing
%        timegrid: equispaced 1D grid (column array) covering all values of 'runtime_over_iter'
% output: t_idx_with_t_on_grid: array with "t_idx_with_t_on_grid(ii)
%                                       = largest index jj with runtime_over_iter(jj) < ii*gridsize + eps"
%         data_with_t_on_grid: rounded data matrix "data( t_idx_with_t_on_grid(ii) )"       


function [t_idx_with_t_on_grid, data_with_t_on_grid] = round_to_timegrid(...
                                                       data, runtime_over_iter, timegrid)
    
    t_idx_with_t_on_grid = zeros( size(timegrid) );

    jj = 1;
    for ii = 1:length(runtime_over_iter) 
        if jj == length(t_idx_with_t_on_grid)
            break
        end
        if runtime_over_iter(ii) < timegrid(jj) + eps
            t_idx_with_t_on_grid(jj) = ii;
        else 
            t_idx_with_t_on_grid(jj+1) = ii;
            jj = jj + 1;
        end
    end

    data_with_t_on_grid = data( t_idx_with_t_on_grid );

end