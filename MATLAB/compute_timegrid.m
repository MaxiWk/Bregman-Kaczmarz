% computes a common timegrid (uniform over all methods and all repeats)
function [num_grid_points, timegrid] = compute_timegrid(runtime_over_iter, num_examples, num_methods)

    % compute gridsize for runtime grid 
    gridsize = 0;
    for example = 1:num_examples
        for method = 1:num_methods
            runtime_over_iter_instance = runtime_over_iter{example, method};
            gridsize = max( [gridsize, runtime_over_iter_instance(1), runtime_over_iter_instance(2:end) - runtime_over_iter_instance(1:end-1)] );
        end
    end
    
    
    % compute common total number of grid points
    num_grid_points = inf;
    for example = 1:num_examples
        for method = 1:num_methods
            runtime_over_iter_instance = runtime_over_iter{example, method};
            num_grid_points = min( num_grid_points, floor(max(runtime_over_iter_instance)/gridsize) );
        end
    end

    timegrid = gridsize * (1:num_grid_points);

end



  
