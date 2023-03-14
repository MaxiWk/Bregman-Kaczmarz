% Algorithm 3
classdef POCS_Simplex_Kaczmarz

    properties
        id = 'POCS';
        %tol_norm_grad = 1e-9 % if norm(grad) < tol_norm_grad, we not perform an update
        % parameters for displaying the error such as color, marker etc
        plot_params = struct('minmaxcolor', [1 0.9 0.9], ...
                                 'quantcolor', [1 0.6 0.6], ...
                                 'linecolor', 'r', ...
                                 'stroke', '--')
    end

    methods

        function vars = update(obj, vars, problem, p)
            i = sampling(p);
            a_i = problem.A(i,:)';
            x = vars.x - (a_i'*vars.x-problem.b(i))/norm(a_i)^2 * a_i;
            vars.x = SimplexProj(x);
        end

    end

end



 