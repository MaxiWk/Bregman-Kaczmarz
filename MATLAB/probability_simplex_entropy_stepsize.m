% solve linesearch problem min_t { beta*t + log(sum(x*exp(-t*alpha)))} 
function t_opt = probability_simplex_entropy_stepsize(x, alpha, beta, solver)

    if solver.make_projection_feasibility_check
        assert( check_if_projection_exists_on_probability_simplex(alpha, beta),...
        'Error in entropy stepsize: Newton equation has no solution in the positive orthant.' );
    end
        
    % first check trivial cases
    if norm(alpha) < 1e-15
        if abs(beta) < 1e-15
            t_opt = 0;
            return;
        else
            error('Error in entropy stepsize: Newton equation has no solution.')
        end
    end
    
    
    % compute stepsize with iterative method
    
    %steps = 0;
    tol_linesearch = solver.tol_linesearch;
    
    t_eucl = ( alpha'*x - beta)/ norm(alpha)^2;   % initial value (t for eucl. proj.)
    t = t_eucl;
    
    exp_term = exp(-t*(alpha));
    sum_x_exp_term = sum(x.*exp_term);
    x_alpha = x.*alpha;
    
    g_1 = beta - sum(x_alpha.*exp_term) / sum_x_exp_term;
    %eval_g_1 = @(t) beta - sum(x_alpha.*exp(-t*alpha)) / sum(x.*exp(-t*alpha));
    %steps = 0;
    
    % debug
    %plot_unit_simplex_entropy_linesearch_function(x, alpha, beta)
    %hold on
    
    switch solver.linesearch_optimizer
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'bisection'
            
            if g_1 < 0
                t_minus = t;
                t_plus = t;
                while g_1 < 0
                    t_plus = t_plus + 1;
                    g_1 = eval_g_1(t_plus);
                end
            else
                t_plus = t;
                t_minus = t;
                while g_1 > 0
                    t_minus = t_minus - 1;
                    g_1 = eval_g_1(t_minus);
                end
            end

            while abs(g_1) > tol_linesearch
                t = 0.5*(t_minus+t_plus);
                g_1 = eval_g_1(t);
                if g_1 > 0
                    t_plus = t;
                else
                    t_minus = t;
                end
                % debug
                %eval_g = @(t) beta*t + log( sum( x.* exp(-t*alpha) ) );
                %plot(t,eval_g(t),'*','color','r')
                %steps = steps + 1;
            end

           %plot(t,eval_g(t),'*','color','r')
           %fprintf('%d steps of bisection method for linesearch \n', steps)
           %fprintf('Stepsize: %f,   initial value (from eucl. proj.): %f \n', t, t_eucl)    
           %eval_g = @(t) beta*t + log( sum( x.* exp(-t*alpha) ) );
           %eval_g_1 = @(t) beta - sum(x_alpha.*exp(-t*alpha)) / sum(x.*exp(-t*alpha));
           %plot(t,eval_g(t),'*','color','green');

         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        case {'Newton', 'regNewton'}
            
            x_exp_term = x.*exp_term;
            alpha_x_exp_term = alpha.*x_exp_term;
            sum_alpha_x_exp_term = sum(alpha_x_exp_term);
            g_1 = beta - sum_alpha_x_exp_term / sum_x_exp_term;
           
        while abs(g_1) > tol_linesearch 
            
            nom_g_2 = sum(alpha.*alpha_x_exp_term) * sum_x_exp_term - sum_alpha_x_exp_term^2;
            denom_g_2 = sum_x_exp_term^2;
            g_2 = nom_g_2 / denom_g_2; 
            
            if strcmp(solver.linesearch_optimizer, 'Newton')
                t = t - g_1/g_2;
            else
                t = t - g_1/(g_2 + solver.Newton_damping_factor * sqrt(abs(g_1)));
            end
            
            exp_term = exp(-t*alpha); %
            x_exp_term = x.*exp_term; %
            alpha_x_exp_term = alpha.*x_exp_term; %
            sum_x_exp_term = sum(x_exp_term); % 
            sum_alpha_x_exp_term = sum(alpha_x_exp_term); %
            g_1 = beta - sum_alpha_x_exp_term / sum_x_exp_term; %
            % debug
            %eval_g = @(t) beta*t + log( sum( x.* exp(-t*alpha) ) );
            %plot(t,eval_g(t),'*','color','red');
            %steps = steps + 1;            
        end

    %plot(t,eval_g(t),'*','color','green');

    % debug
    %fprintf('%d steps of Newton`s method for linesearch \n', steps)
    %fprintf('Stepsize: %f,   initial value (from eucl. proj.): %f \n', t, t_eucl)

    %if steps > 30
    %    fprintf('More than 30 steps of Newton`s method for linesearch (%d) \n', steps)
    %end

    end
    
    
    t_opt = t;
   

    
end