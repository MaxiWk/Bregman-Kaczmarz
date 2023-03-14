% Nonlinear Kaczmarz method based on euclidean projections 
classdef NK  

    properties
        id = 'Nonlinear Kaczmarz';
        tol_norm_grad = 1e-9 % if norm(grad) < tol_norm_grad, we not perform an update
        % parameters for displaying the error such as color, marker etc
        plot_params = struct('minmaxcolor', [1 0.9 0.9], ...
                                 'quantcolor', [1 0.6 0.6], ...
                                 'linecolor', 'r', ...
                                 'stroke', '--')
    end

    methods

        function vars = update(obj, vars, problem, p)
            i = sampling(p);
            compute_gradient = true;
            [F_i_x, grad_F_i_x] = problem.eval_f_and_grad_f(i, vars.x, compute_gradient);
            norm_grad_F_i_x_sqr = norm(grad_F_i_x)^2;
            if norm(norm_grad_F_i_x_sqr) > obj.tol_norm_grad
                vars.x = vars.x - F_i_x/norm_grad_F_i_x_sqr *grad_F_i_x;
            end

        end

    end

end


