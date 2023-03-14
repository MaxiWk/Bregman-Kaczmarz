classdef sparse_solutions_of_quadratic_equations    
    
    properties 
        system_size % struct(d,n) with d: number of dimensions, n: number of equations
        A % matrices Ai in the system 
        b % vectors bi in the system 
        c % scalars ci in the system
        sol % particular exact solution of the system
        init_vars % initial values for iterates of the solvers
    end

    
    methods


        function obj = sparse_solutions_of_quadratic_equations(problem_data)

                obj.system_size = problem_data.system_size;
                
                rand('state', problem_data.random_seed_for_problem.current);
                randn('state', problem_data.random_seed_for_problem.current);
                
                obj.sol = obj.setup_sol(problem_data.system_size.d, problem_data.sol_properties); 
                [obj.A, obj.b] = obj.setup_Ab(problem_data.system_size, problem_data.distr_params);
                
                % choose c from A,b and sol such that f(sol) = 0
                c = zeros(problem_data.system_size.n,1);
                for eq = 1:problem_data.system_size.n
                    c(eq) = -obj.sol'* (0.5*obj.A(:,:,eq)*obj.sol + obj.b(:,eq));
                end 
                obj.c = c;
                
                if isfield(problem_data, 'random_seed_for_init')
                    rand('state', problem_data.random_seed_for_init.current);
                    randn('state', problem_data.random_seed_for_init.current);
                end
                obj.init_vars = obj.setup_init_vars(problem_data.system_size.d, problem_data.lambda, ...
                                                problem_data.init_properties);

        end




       function varargout = eval_f_and_grad_f(obj, i, x, compute_gradient)

            Ai_x = obj.A(:,:,i) * x;
            bi = obj.b(:,i);
            F_i_x = x' * (0.5*Ai_x + bi) + obj.c(i);
            varargout{1} = F_i_x;

            if compute_gradient
                grad_F_i_x = Ai_x + bi; 
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
            end
        end




        function sol = setup_sol(obj, d, sol_properties)

            assert(isfield(sol_properties, 'type'))
            assert(isfield(sol_properties, 'nnz'))
        
            switch sol_properties.type
                case 'sparse_randn'
                    sol_nonzero_part = randn(sol_properties.nnz, 1);
                case 'sparse_rand'
                    assert(isfield(sol_properties, 'min_nonzero'))
                    assert(isfield(sol_properties, 'max_nonzero'))
                    width = sol_properties.max_nonzero - sol_properties.min_nonzero;
                    assert(with > 0, 'Choose sol_properties.min_nonzero < sol_properties.max_nonzero!')
                    sol_nonzero_part = sol_properties.min_nonzero...
                                            + width * rand(sol_properties.nnz, 1);
                otherwise 
                    disp('This sol_properties.type is not implemented!')
            end
        
            sol = zeros(d, 1);
            rand_perm = randperm(d);
            rand_pos = rand_perm(1:sol_properties.nnz);   
            sol(rand_pos) = sol_nonzero_part;
    
        end




        function init_vars = setup_init_vars(obj, d, lambda, init_properties)

            switch init_properties.type
                case 'dual_randn'
                    init_vars.xstar = randn(d, 1);
                    init_vars.x = grad_phistar(init_vars.xstar, lambda);
                case 'zero'
                    init_vars.xstar = zeros(d,1);
                    init_vars.x = grad_phistar(init_vars.xstar, lambda);
                case 'local'
                    assert(isfield(init_properties, 'init_dist_x'))
                    rand_direction = rand(d, 1);
                    rand_direction = rand_direction/norm(rand_direction);
                    init_vars.x = obj.sol ...
                                        + init_properties.init_dist_x * rand_direction;
                    init_vars.xstar = grad_phistar_inv(init_vars.x, lambda);
                otherwise
                    disp('This init_type is not implemented!')
            end

        end





        function [A,b] = setup_Ab(obj, system_size, distr_params)

            d = system_size.d;
            n = system_size.n;
            A_distr = distr_params.A_distr;
            b_distr = distr_params.b_distr;
        
            A = zeros(d, d, n);
            b = zeros(d, n);
        
            for eq = 1:n
        
                switch A_distr
                    case 'randn'
                        A_for_single_eq = randn(d, d);
                    case 'randn_plus_diag'
                        assert(isfield(distr_params, 'reg_sing_val'))
                        A_for_single_eq = randn(d, d) + distr_params.reg_sing_val * eye(d, d);
                    case 'randn_symm'
                        A_for_single_eq = randn(d, d);
                        A_for_single_eq = 0.5*(A_for_single_eq + A_for_single_eq'); % symmetrize A
                    case 'rand'
                        A_for_single_eq = rand(d, d);
                    otherwise 
                        disp('This A_distr is not implemented!');
                end
        
                switch b_distr
                    case 'randn'
                        b_for_single_eq = randn(d, 1);
                    case 'rand'
                        b_for_single_eq = rand(d, 1);
                    otherwise
                        disp('This b_distr is not implemented!')
                end
        
                A(:,:,eq) = A_for_single_eq;
                b(:,eq) = b_for_single_eq;
        
            end
            
        end



    end
end