% Algorithm 2 for the case of sparse solutions of quadratic equations 
classdef rNBK  % relaxed NBK method with grad_phistar = soft shrinkage 

    properties
        id = 'Relaxed Sparse NBK'
        lambda
        tol_norm_grad = 1e-9 % if norm(grad) < tol_norm_grad, we not perform an update
        % parameters for displaying the error such as color, marker etc
        plot_params = struct('minmaxcolor', [0.9 0.9 1], ...
                                 'quantcolor', [0.6 0.6 1], ...
                                 'linecolor', 'b', ...
                                 'stroke', ':')
    end

    methods

        function obj = rNBK(lambda)
            obj.lambda = lambda;
        end

        function vars = update(obj, vars, problem, p)
            i = sampling(p);
            compute_gradient = true;
            [F_i_x, grad_F_i_x] = problem.eval_f_and_grad_f(i, vars.x, compute_gradient);
            norm_grad_F_i_x_sqr = norm(grad_F_i_x)^2;
            if norm(norm_grad_F_i_x_sqr) > obj.tol_norm_grad
                vars.xstar = vars.xstar - F_i_x/norm_grad_F_i_x_sqr *grad_F_i_x;
                vars.x = grad_phistar(vars.xstar, obj.lambda);
            end

        end

    end

end


