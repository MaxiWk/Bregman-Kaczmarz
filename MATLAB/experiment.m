% specifies experiment and how to run it
classdef experiment

    properties
        problem_class % problem which is to be solved (e.g. sparse_solutions_of_quadratic_equations)
        problem_data % struct which defines the problem (constructor takes these parameters and outputs new random problem)
        solvers % cell array of solvers (NBK, rNBK, ...)    
        num_solvers % length(solvers)
        num_examples % number of random instances       
        stopping % stopping class
        result_types % types of results which shall be reported, e.g. 'residual', 'nnz'
        p % sampling probabilities for the equations
        random_seeds % struct containing seeds for random generator, 
                     % seperately for problem generation, initial values
                     % and the algorithm (all three require random numbers)
                     % ('start': an initial seed, 'increment' decides if 
                     % different random numbers shall appear in each 
                     % random experiment, 'current' is the current seed)
    end
    

    methods


        function obj = experiment(problem_class, problem_data, solvers, num_examples, stopping, result_types, p, random_seeds)

            arguments 
                problem_class 
                problem_data 
                solvers cell
                num_examples int32
                stopping 
                result_types cell
                p (1,:) double
                random_seeds struct
            end

            obj.problem_class = problem_class;
            obj.problem_data = problem_data;
            obj.solvers = solvers; 
            obj.num_solvers = length(solvers);
            obj.num_examples = num_examples;
            obj.stopping = stopping;
            obj.result_types = result_types;
            obj.p = p;

            if ~isfield(random_seeds, 'random_seed_for_problem')
                random_seeds.random_seed_for_problem = struct('start', 0, 'increment', 1);
            end

            if ~isfield(random_seeds, 'random_seed_for_solver')
                random_seeds.random_seed_for_solver = struct('start', 0, 'increment', 1);
            end            

            obj.random_seeds = random_seeds;
            
        end



        function stats = run(obj)

            % initialization of arrays which store the errors
            num_saved_iters = floor(obj.stopping.maxiter / obj.stopping.iter_save);
            all_results = init_all_results(obj.result_types, obj.stopping, ...
                                            num_saved_iters, ...
                                            obj.num_examples, obj.num_solvers);

            % run solvers 
            for example_counter = 1:obj.num_examples

                    fprintf('Example %d/%d \n', example_counter, obj.num_examples)

                for solver_counter = 1:obj.num_solvers

                    % set up problem                    
                    obj.random_seeds.random_seed_for_problem...
                            = set_seed(obj.random_seeds.random_seed_for_problem, example_counter); 
                    obj.problem_data.random_seed_for_problem = obj.random_seeds.random_seed_for_problem;
                    if isfield(obj.random_seeds, 'random_seed_for_init')
                        obj.random_seeds.random_seed_for_init...
                            = set_seed(obj.random_seeds.random_seed_for_init, example_counter);   
                        obj.problem_data.random_seed_for_init = obj.random_seeds.random_seed_for_init;
                    end
                    problem = obj.problem_class(obj.problem_data); 

                    % choose solver
                    solver = obj.solvers{solver_counter};
                    if isfield(obj.random_seeds, 'random_seed_for_solvers')
                        obj.random_seeds.random_seed_for_solvers = set_seed(obj.random_seeds.random_seed_for_solvers, example_counter);
                        rand('state', obj.random_seeds.random_seed_for_solvers.current)
                    end
                    fprintf('Running %s', solver.id)

                    % run solver 
                    tic
                    new_results = run_solver(problem, solver, obj.p, obj.stopping, obj.result_types);
                    time = toc;

                    % report result
                    fprintf('\t Computation time: %f seconds \n', time)
                    if any(strcmp(obj.result_types, 'total_runtime'))
                        new_results.total_runtime = time;
                    end
                    all_results = report_new_result(all_results, ...
                                      new_results, obj.stopping, example_counter, solver_counter);

                end
            end 

            % compute error statistics 
            stats = compute_stats_for_all_results(... 
                            all_results, obj.stopping, obj.num_solvers, obj.num_examples);
 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % function: update random seed
            function seed_struct = set_seed(seed_struct, num_current_example)
                if num_current_example == 1
                    seed_struct.current = seed_struct.start;
                else
                    seed_struct.current = seed_struct.current + seed_struct.increment;
                end
            end
            
            % function: write the errors of the method into the
            % corresponding matrices in full_error_data 
            function all_results = report_new_result(all_results, results,...
                                            stopping, example_counter, solver_counter)
                fn = fieldnames(all_results);
                for field_num = 1:length(fn)
                    key = fn{field_num};
                    val = all_results.(key);
                    if any(strcmp(key, {'nnz', 'total_runtime'})) % error per run
                       val(example_counter, solver_counter) = results.(key);
                    elseif any(strcmp(key, {'residual', 'dist_to_sol', 'runtime_over_iter'})) % error per run and saved iteration
                       if isa(stopping, 'iteration_stopping')
                          val(:, example_counter, solver_counter) = results.(key);  
                       elseif isa(stopping, 'runtime_stopping')
                          val{example_counter, solver_counter} = results.(key); 
                       else 
                          disp('This stopping is not implemented!')
                       end
                    end
                    all_results.(key) = val;
                end
            end

        end

    end

end