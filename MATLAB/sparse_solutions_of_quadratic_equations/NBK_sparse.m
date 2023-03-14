% Algorithm 1 for the case of sparse solutions of quadratic equations 
classdef NBK_sparse  % NBK with DGF phi(x) = lambda*||x||_1 + 1/2*||x||_2^2 

    properties
        id = 'Sparse NBK';
        lambda % parameter in DGF
        tol_norm_grad = 1e-9 % if norm(grad) < tol_norm_grad, we not perform an update
        % parameters for displaying the error such as color, marker etc
        plot_params = struct('minmaxcolor', [0.9 1 0.9], ...
                                 'quantcolor', [0.6 1 0.6], ...
                                 'linecolor', 'g', ...
                                 'stroke', '-')
    end
    
    methods

        function obj = NBK_sparse(lambda)
            obj.lambda = lambda;
        end

        function vars = update(obj, vars, problem, p)
            i = sampling(p);
            compute_gradient = true;
            [F_i_x, grad_F_i_x] = problem.eval_f_and_grad_f(i, vars.x, compute_gradient);
            if norm(grad_F_i_x) > obj.tol_norm_grad
                beta = grad_F_i_x'* vars.x - F_i_x;                
                [vars, ~] = linesearch_shrinkage(vars, grad_F_i_x, beta, obj.lambda);
            end

        end

    end

end


