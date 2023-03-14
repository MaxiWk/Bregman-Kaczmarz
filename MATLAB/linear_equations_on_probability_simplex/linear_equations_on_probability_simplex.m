classdef linear_equations_on_probability_simplex    
    
    properties 
        system_size % struct(d,n) with d x n: size of system matrix A
        A % matrix A in the system 
        b % vector b in the system 
        sol % particular exact solution of the system
        init_vars % initial values for iterates of the solvers
    end

    
    methods


        function obj = linear_equations_on_probability_simplex(problem_data)

                obj.system_size = problem_data.system_size;

                rand('state', problem_data.random_seed_for_problem.current);
                randn('state', problem_data.random_seed_for_problem.current);
                obj.sol = obj.setup_sol(problem_data.system_size.d, problem_data.sol_distr);
                obj.A = obj.setup_A(problem_data.system_size, problem_data.A_distr);

                obj.b = obj.A * obj.sol; 
                if isfield(problem_data, 'random_seed_for_init')
                    rand('state', problem_data.random_seed_for_init.current);
                    randn('state', problem_data.random_seed_for_init.current);
                end
                obj.init_vars = obj.setup_init_vars(problem_data.system_size.d, problem_data.init_properties);
                                                
        end



        function varargout = eval_f_and_grad_f(obj, i, x, compute_gradient)

            F_i_x = obj.A(i,:) * x - obj.b(i);
            varargout{1} = F_i_x;

            if compute_gradient
                grad_F_i_x = obj.A(i,:)';
                varargout{2} = grad_F_i_x;                
            end

        end



        function error = compute_error(obj, vars, result_type)
            switch result_type
                case 'dist_to_sol'
                    error = norm(vars.x - obj.sol);
                case 'residual'
                    res = 0;
                    compute_gradient = false;                         
                    for i = 1:obj.system_size.n
                        F_i_x = obj.eval_f_and_grad_f(i, vars.x, compute_gradient);
                        res = res + F_i_x^2;
                    end                        
                    error = sqrt(res);  
                    error = error / norm(obj.b);
            end
        end



        function sol = setup_sol(obj, d, sol_distr)

            switch sol_distr
                case 'random_uniform'
                    sol = random_point_on_probability_simplex(d);
                otherwise 
                    disp('This sol_properties.type is not implemented!')
            end
        
        end






        function init_vars = setup_init_vars(obj, d, init_properties)
        
            switch init_properties.type
                case 'random_uniform'
                    init_vars.x = random_point_on_probability_simplex(d);
                case 'central'
                    init_vars.x = ones(d,1)/d;
                otherwise
                    disp('This init_type is not implemented!')
            end
        
        end






         function A = setup_A(obj, system_size, A_distr)
        
            d = system_size.d;
            n = system_size.n;
        
            switch A_distr.type
                case 'randn'
                    A = randn(n, d);
                case 'uniform'
                    width = A_distr.max_A - A_distr.min_A;
                    A = A_distr.min_A + rand(n, d) * width;
                case 'split_uniform'
                    rand_sign = 2 * randi([0 1], n, d) - 1;
                    A = rand_sign .* ( A_distr.abs_min_A + ...
                    (A_distr.abs_max_A - A_distr.abs_min_A) * rand(n, d) );
                otherwise 
                    disp('This A_distr is not implemented!');
            end
        
         end

    end

end