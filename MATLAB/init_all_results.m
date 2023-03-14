% generates empty arrays which store error quantities, which are specified
% in cell array 'result_types'
function all_results = init_all_results(result_types, stopping, ...
                                             num_result_to_save_per_solver_and_example, ...
                                             num_examples, num_solvers)

    assert(iscell(result_types));

    all_results = struct();

    % error quantities whose size do not depend on the number of iterations
    for ii = 1:length(result_types)
        if any(strcmp(result_types{ii}, {'nnz', 'total_runtime'})) 
           all_results.(result_types{ii}) = zeros(num_examples, num_solvers);
        end
    end

    % error quantities whose size depend on the number of iterations, here:
    % stopping after a specified number of iterations
    % (then, since errors are saved after a specific number of iterations,
    %  the number of stored values is independent on the random example)
    if isa(stopping, 'iteration_stopping')
        for ii = 1:length(result_types)
            if any(strcmp(result_types{ii}, {'residual', 'dist_to_sol'}))  
                all_results.(result_types{ii}) = zeros(num_result_to_save_per_solver_and_example, ...
                                             num_examples, num_solvers);
            end
        end


    % error quantities whose size depend on the number of iterations, here:
    % stopping after a specified computation time
    % (here, the number of stored values will depend on the random example)
    elseif isa(stopping, 'runtime_stopping')
           all_results.runtime_over_iter = cell(num_examples, num_solvers);                   
           for ii = 1:length(result_types)
               if any(strcmp(result_types{ii}, {'residual', 'dist_to_sol'}))
                  all_results.(result_types{ii}) = cell(num_examples, num_solvers);
               end   
           end

    else
        disp('This stopping type is not implemented!')
    end

end