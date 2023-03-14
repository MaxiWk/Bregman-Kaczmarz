% Algorithm 1 for the case of linear equations on the probability simplex
classdef NBK_Simplex_Kaczmarz  % NBK with negative DGF on simplex 

    properties
        id = 'NBK';
        %tol_norm_grad = 1e-9 % if norm(grad) < tol_norm_grad, we not perform an update
        % parameters for displaying the error such as color, marker etc
        make_projection_feasibility_check
        num_empty_intersections
        tol_linesearch
        linesearch_optimizer
        Newton_damping_factor
        plot_params = struct('minmaxcolor', [0.9 1 0.9], ...
                             'quantcolor', [0.6 1 0.6], ...
                             'linecolor', 'g', ...
                             'stroke', '-')
    end 
    
    methods

        function obj = NBK_Simplex_Kaczmarz(init_struct)
            obj.make_projection_feasibility_check = init_struct.make_projection_feasibility_check;
            obj.tol_linesearch = init_struct.tol_linesearch;
            obj.linesearch_optimizer = init_struct.linesearch_optimizer;
            if strcmp(obj.linesearch_optimizer, 'regNewton')
                assert(isfield(init_struct, 'Newton_damping_factor'))
                obj.Newton_damping_factor = init_struct.Newton_damping_factor;
            end
        end

        function vars = update(obj, vars, problem, p)
            i = sampling(p);
            a_i = problem.A(i,:)';
            t = probability_simplex_entropy_stepsize(vars.x, a_i, problem.b(i), obj);
            x = vars.x .* exp(-t * a_i);
            vars.x = x / sum(x);
        end

    end 

end


