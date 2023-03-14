% Algorithm 2 for the left stochastic decomposition problem 
classdef LSD_rNBK

    properties
        id = 'rNBK'; % 'projected nonlinear Kaczmarz'
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
            j = sampling(p);
        
            x_1 = vars.P(:,i);
        
            if i == j
                grad_1 = 2 * x_1;
                F_i_X = x_1' * x_1 - problem.A(i,i);
                t = F_i_X / norm(grad_1,'inf')^2;
                x_1 = x_1 .* exp(-t * grad_1); 
                vars.P(:,i) = x_1 / sum(x_1);  
            else
                x_2 = vars.P(:,j);
                grad_1 = x_2;
                grad_2 = x_1;
                F_i_X = x_1' * x_2 - problem.A(i,j);
                t = F_i_X / (norm(grad_1,'inf')^2 + norm(grad_2,'inf')^2);
                x_1 = x_1 .* exp(-t * grad_1);
                vars.P(:,i) = x_1 / sum(x_1);  
                x_2 = x_2 .* exp(-t * grad_2);
                vars.P(:,j) = x_2 / sum(x_2);  
                
            end

        end

    end

end



 