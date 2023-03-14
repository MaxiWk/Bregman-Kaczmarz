% Algorithm 4 ('projected nonlinear Kaczmarz method') 
classdef LSD_PNK    

    properties
        id = 'PNK'; 
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
            j = sampling(p);
        
            x_1 = vars.P(:,i);
        
            if i == j
                grad_1 = 2 * x_1;
                F_i_X = x_1' * x_1 - problem.A(i,i);
                t = F_i_X / norm(grad_1)^2;
                x_1 = x_1 - t * grad_1;
                vars.P(:,i) = SimplexProj(x_1);        
            else
                x_2 = vars.P(:,j);
                grad_1 = x_2;
                grad_2 = x_1;
                F_i_X = x_1' * x_2 - problem.A(i,j);
                t = F_i_X / (norm(grad_1)^2 + norm(grad_2)^2);
                x_1 = x_1 - t * grad_1;
                vars.P(:,i) = SimplexProj(x_1);
                x_2 = x_2 - t * grad_2;
                vars.P(:,j) = SimplexProj(x_2);
                
            end

        end

    end

end



 