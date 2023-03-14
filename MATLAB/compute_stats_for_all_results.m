% computes min, max, median and 25%- and 75% quantile w.r.t. all random instances 
% output: struct 'stats', where 
%         - if key == 'nnz, total_runtime' (errors per run),
%               val is a cell array containing 5 arrays of size num_solvers
%         - if key == 'res, dist_to_sol, ...' (errors per run and saved iterations),
%               val is a cell array containing 5 matrices of size num_saved_iters x num_solvers
%         - in case of runtime_stopping, an entry 'timegrid' is added for
%               plotting etc (here we round all errors to a common
%               timegrid for averaging over different random instances)
function stats = compute_stats_for_all_results(all_results, stopping, num_solvers, num_examples)


    % here, stats are computed for errors which are given per run (nnz,
    % total_runtime, ...)
    fn = fieldnames(all_results);
    for field_num = 1:length(fn)
        key = fn{field_num};
        if any(strcmp(key, {'nnz', 'total_runtime'})) 
           val = all_results.(key);
           if size(val,1) == 1 % if only one example, nothing to compute
              [stats.(key).mins, stats.(key).maxs, stats.(key).medians, ...
               stats.(key).quant25s, stats.(key).quant75s] = deal(val);
               
           else
               stats.(key).mins = min(val,[],1);
               stats.(key).maxs = max(val,[],1);
               stats.(key).medians = median(val,1);
               stats.(key).quant25s = quantile(val,0.25,1);
               stats.(key).quant75s = quantile(val,0.75,1);
           end
        end
    end


    % preprocessing in case of runtime_stopping: here, we setup a 
    % common time grid 'rounded_runtimes' and round all error quantities
    % correspondingly
    if isa(stopping, 'runtime_stopping') 
       [num_grid_points, stats.timegrid] = compute_timegrid(...
           all_results.runtime_over_iter, num_examples, num_solvers);
       fn = fieldnames(all_results);
       for field_num = 1:length(fn)
           key = fn{field_num};           
           val = all_results.(key);
           if any(strcmp(key, {'residual', 'dist_to_sol', 'velocity'}))
               rounded_val = zeros(num_grid_points, num_examples, num_solvers);
               for example_counter = 1:num_examples
                   for solver_counter = 1:num_solvers
                       [~, rounded_val(:, example_counter, solver_counter)] = round_to_timegrid(val{example_counter, solver_counter}, ...
                                            all_results.runtime_over_iter{example_counter, solver_counter}, ...
                                            stats.timegrid);     
                   end 
               end
               all_results.(key) = rounded_val;
           end
       end
    end


   % here, stats are computed for errors which are given per saved iteration (res,
   % dist_to_sol, ...)
    fn = fieldnames(all_results);
    for field_num = 1:length(fn)
       key = fn{field_num};
       if any(strcmp(key, {'residual', 'dist_to_sol'})) 
           val = all_results.(key);
           if size(val,2) == 1   % if only one example, nothing to compute 
               [stats.(key).mins, stats.(key).maxs, stats.(key).medians, ...
                stats.(key).quant25s, stats.(key).quant75s]...
                    = deal(reshape(val, size(val,1), size(val,3)));
           else
              for solver_counter = 1:size(val,3) % loop over solvers
                  stats.(key).mins(:,solver_counter) = min(val(:,:,solver_counter), [], 2);
                  stats.(key).maxs(:,solver_counter) = max(val(:,:,solver_counter), [], 2);
                  stats.(key).medians(:,solver_counter) = median(val(:,:,solver_counter),2);
                  stats.(key).quant25s(:,solver_counter) = quantile(val(:,:,solver_counter),0.25,2); 
                  stats.(key).quant75s(:,solver_counter) = quantile(val(:,:,solver_counter),0.75,2);
              end
           end  
       end
    end

end
