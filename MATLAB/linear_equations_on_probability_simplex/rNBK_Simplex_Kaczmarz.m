% Algorithm 2 for the case of linear equations on the probability simplex
classdef rNBK_Simplex_Kaczmarz

    properties
        id = 'rNBK';
        %tol_norm_grad = 1e-9 % if norm(grad) < tol_norm_grad, we not perform an update
        % parameters for displaying the error such as color, marker etc
        plot_params = struct('minmaxcolor', [0.9 0.9 1], ...
                                 'quantcolor', [0.6 0.6 1], ...
                                 'linecolor', 'b', ...
                                 'stroke', ':')   
    end

    methods
        function vars = update(obj, vars, problem, p)
            i = sampling(p);
            a_i = problem.A(i,:)';
            b_i = problem.b(i);
            t = (a_i'*vars.x-b_i)/norm(a_i, 'Inf')^2;
            x = vars.x .* exp(-t * a_i);
            vars.x = x / sum(x);           
        end
    end

end