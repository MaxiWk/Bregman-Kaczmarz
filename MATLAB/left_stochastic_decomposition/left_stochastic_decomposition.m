classdef left_stochastic_decomposition    
    
    properties 
        system_size % struct(r,m) with d: number of rows, n: number of columns 
        A % matrix A in the system 
        sol % particular exact solution of the system
        init_vars % initial values for iterates of the solvers
    end

    
    methods


        function obj = left_stochastic_decomposition(problem_data)

            obj.system_size = problem_data.system_size; 

            rand('state', problem_data.random_seed_for_problem.current);
            randn('state', problem_data.random_seed_for_problem.current);          
            obj.sol = obj.setup_sol(problem_data.system_size, problem_data.sol_distr);
            obj.A = obj.sol' * obj.sol;

            if isfield(problem_data, 'random_seed_for_init')
                rand('state', problem_data.random_seed_for_init.current);
                randn('state', problem_data.random_seed_for_init.current);
            end
            obj.init_vars = obj.setup_init_vars(problem_data.system_size, problem_data.init_properties);
                                                
        end



        function varargout = eval_f_and_grad_f(obj, i, j, P, compute_gradient)

            F_i_x = P(:,i)'*P(:,j) - obj.A(i,j);
            varargout{1} = F_i_x;

            if compute_gradient
               grad_F_i_X = ((1:pars.m)==i).*P(:,j) + ((1:pars.m)==j).*P(:,i);
               varargout{2} = grad_F_i_X;             
            end

        end



        function error = compute_error(obj, vars, result_type)
            switch result_type 
                case 'dist_to_sol'
                    error = norm(vars.P - obj.sol, 'fro');
                case 'residual'
                    res = 0;
                    compute_gradient = false;                         
                    for i = 1:obj.system_size.m
                        for j = 1:obj.system_size.m
                            F_i_x = obj.eval_f_and_grad_f(i, j, vars.P, compute_gradient);
                            res = res + F_i_x^2;
                        end
                    end                        
                    error = sqrt(res);  
            end
        end



        function sol = setup_sol(obj, system_size, sol_distr)

            sol = zeros(system_size.r, system_size.m);

            switch sol_distr.type 
                case 'random_uniform'
                    for i = 1:system_size.m
                        sol(:,i) = random_point_on_probability_simplex(system_size.r);
                    end
                case 'random_redundant'
                    assert(isfield(sol_distr, 'redundancy'))
                    assert(isfield(sol_distr, 'redundancy_noiselev'))
                    assert(sol_distr.redundancy >= 0 && sol_distr.redundancy <= 1);                        
                    pivot = round((1-sol_distr.redundancy) * system_size.m);
                    % initialize the first columns independently 
                    for i = 1:pivot
                        sol(:, i) = random_point_on_probability_simplex(pars.r); 
                    end
                    % create noisy copies of same block size until end of array
                    last_filled_col = pivot;
                    while last_filled_col < system_size.m
                        num_cols_to_fill = min(pivot, system_size.m - last_filled_col);
                        left = last_filled_col + 1;
                        right = last_filled_col + num_cols_to_fill;
                        sol(:, left:right) = sol(:, 1:num_cols_to_fill)...                   
                                    + sol_distr.redundancy_noiselev ...
                                    * rand(system_size.r, num_cols_to_fill);
                        last_filled_col = last_filled_col + num_cols_to_fill;
                    end
                    % project the columns of sol back to the simplex
                    sol = max(sol, 0);
                    sol = bsxfun(@rdivide, sol, sum(sol,1));                    
                otherwise 
                    disp('This sol_properties.type is not implemented!')
            end
        
        end






        function init_vars = setup_init_vars(obj, system_size, init_properties)
        
            init_vars.num_empty_intersections = 0;
            switch init_properties.type
                case 'random_uniform'
                    init_vars.P = zeros(system_size.r, system_size.m);
                    for i = 1:system_size.m
                        init_vars.P(:,i) = random_point_on_probability_simplex(system_size.r);
                    end
                case 'central'
                    init_vars.P = ones(system_size.r, system_size.m)/system_size.r;
                otherwise
                    disp('This init_type is not implemented!')
            end
        
        end




    end

end