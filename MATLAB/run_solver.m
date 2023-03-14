% runs 'solver' on 'problem' with 'stopping' condition 
% equations are sampled randomly according to the probability vector 'p'
% errors specified as strings in 'result_types' are reported in 'results'
function results = run_solver(problem, solver, p, stopping, result_types)

    stop = false;
    iter = 0; 
    iter_plot = 1;
    
    vars = problem.init_vars;
    
    while ~stop

        iter = iter + 1;

        vars = solver.update(vars, problem, p);        

       % check stopping condition
       if isa(stopping, 'iteration_stopping')
           stop = (iter > stopping.maxiter);
       elseif isa(stopping, 'runtime_stopping')
           stop = (toc > stopping.max_time);
       elseif isa(stopping, 'residual_stopping')
           res = problem.compute_errors(vars, 'res');
           stop = (res > stopping.max_res);
       else
           disp('This stopping is not implemented!') 
       end



        % each 'iter_save' many iterations, report errors 
        if mod(iter, stopping.iter_save) == 0 

            if iter_plot > round(stopping.maxiter / stopping.iter_save)
                continue
            end

            for ii = 1:length(result_types)
                result_type = result_types{ii};
                if strcmp(result_type, 'nnz')  
                    continue % will be reported once at the end
                end
                if strcmp(result_type, 'num_empty_intersections') 
                    continue % will be reported once at the end
                end                        
                if strcmp(result_type, 'total_runtime') 
                    continue % will be reported once at the end in 'experiment.m'
                end 
                if isa(stopping, 'residual_stopping') && strcmp(result_type, 'residual') 
                    continue % already computed
                end        
                if isa(stopping, 'runtime_stopping')
                    results.runtime_over_iter(iter_plot) = toc;
                end
                results.(result_type)(iter_plot) = problem.compute_error(vars, result_type); 
            end
            iter_plot = iter_plot + 1;

        end

    end 


    % report sparsity of last iterate 
    if any(strcmp(result_types, 'nnz'))
       tol_sparsity = 1e-5;
       results.nnz = length(find((abs(vars.x) > tol_sparsity)));
    end

    % report number of iterations in which the intersection condition was
    % not fulfilled
    if any(strcmp(result_types, 'num_empty_intersections'))
       results.num_empty_intersections = vars.num_empty_intersections;
    end

end

