% solve linesearch problem min_t { beta*t + log(sum(x_1*exp(-t*alpha_1))) + log(sum(x_2*exp(-t*alpha_2))) } 
function t_opt = product_probability_simplex_entropy_stepsize(x_1, x_2, alpha_1, alpha_2, beta, solver)

    %plot_double_unit_simplex_entropy_linesearch_function(x_1, x_2, alpha_1, alpha_2, beta)

    if solver.make_projection_feasibility_check
        assert( check_if_projection_exists_on_product_probability_simplex(alpha_1, alpha_2, beta),...
        'Error in entropy stepsize: Newton equation has no solution on the constraint set.' );
    end 
        
    % first check trivial cases
    if sqrt(norm(alpha_1)^2+norm(alpha_2)^2) < 1e-15 
        if abs(beta) < 1e-15
            t_opt = 0;
            return;
        else
            error('Error in entropy stepsize: Newton equation has no solution.')
        end
    end
    
      
    tol_linesearch = solver.tol_linesearch;
    
    t_eucl = ( alpha_1'*x_1 + alpha_2'*x_2 - beta)/ (norm(alpha_1)^2 + norm(alpha_2)^2);   % initial value (t for eucl. proj.)
    t = t_eucl;

    % first step
    if solver.stabilize_exp % more expensive operations, but more stable expressions (if needed)

        alpha_1_max = max(alpha_1);
        alpha_2_max = max(alpha_2);

        exp_term_1 = exp(t*(alpha_1_max-alpha_1));
        exp_term_2 = exp(t*(alpha_2_max-alpha_2));
        sum_x_exp_term_1 = sum(x_1.*exp_term_1);
        sum_x_exp_term_2 = sum(x_2.*exp_term_2);
        x_alpha_1 = x_1.*alpha_1;
        x_alpha_2 = x_2.*alpha_2;
        
        g_1 = beta - sum(x_alpha_1.*exp_term_1) / sum_x_exp_term_1 - sum(x_alpha_2.*exp_term_2) / sum_x_exp_term_2;
    
        %eval_g = @(t) beta*t + log( sum( x_1.* exp(-t*alpha_1) ) )...
        %                     + log( sum( x_2.* exp(-t*alpha_2) ) );    
        eval_g_1 = @(t) beta - sum(x_alpha_1.*exp(-t*alpha_1)) / sum(x_1.*exp(-t*alpha_1))...
                         - sum(x_alpha_2.*exp(-t*alpha_2)) / sum(x_2.*exp(-t*alpha_2));

    else 

        exp_term_1 = exp(-t*alpha_1);
        exp_term_2 = exp(-t*alpha_2);
        sum_x_exp_term_1 = sum(x_1.*exp_term_1);
        sum_x_exp_term_2 = sum(x_2.*exp_term_2);
        x_alpha_1 = x_1.*alpha_1;
        x_alpha_2 = x_2.*alpha_2;
        
        g_1 = beta - sum(x_alpha_1.*exp_term_1) / sum_x_exp_term_1 - sum(x_alpha_2.*exp_term_2) / sum_x_exp_term_2;
    
        %eval_g = @(t) beta*t + log( sum( x_1.* exp(-t*alpha_1) ) )...
        %                     + log( sum( x_2.* exp(-t*alpha_2) ) );    
        eval_g_1 = @(t) beta - sum(x_alpha_1.*exp(-t*alpha_1)) / sum(x_1.*exp(-t*alpha_1))...
                         - sum(x_alpha_2.*exp(-t*alpha_2)) / sum(x_2.*exp(-t*alpha_2));

    end

    %steps = 0;
    
    % debug
    % plot_double_unit_simplex_entropy_linesearch_function(x_1, x_2, alpha_1, alpha_2, beta)
    % hold on
    
    switch solver.linesearch_optimizer
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'bisection'
            
            % find initial boundary points
            if g_1 < 0
                t_minus = t;
                t_plus = t;
                while g_1 < 0
                    t_plus = t_plus + 1;
                    g_1 = eval_g_1(t_plus);
                    if isprop(solver, 'max_abs_t')
                        assert( abs(t_plus) < solver.max_abs_t )
                    end 
                end
            else
                t_plus = t;
                t_minus = t;
                while g_1 > 0
                    t_minus = t_minus - 1;
                    g_1 = eval_g_1(t_minus);
                    if isfield(pars, 'max_abs_t')
                        assert( abs(t_minus) < solver.max_abs_t )
                    end                       
                end
            end        

            % iterations
            while abs(g_1) > tol_linesearch
                t = 0.5*(t_minus+t_plus);
                g_1 = eval_g_1(t);
                if g_1 > 0
                    t_plus = t;
                else
                    t_minus = t;
                end
                % debug
                %plot(t,eval_g(t),'*','color','r')
                %steps = steps + 1;
            end

           %plot(t,eval_g(t),'*','color','r')
           %fprintf('%d steps of bisection method for linesearch \n', steps)
           %fprintf('Stepsize: %f,   initial value (from eucl. proj.): %f \n', t, t_eucl)    
           %plot(t,eval_g(t),'*','color','green');

         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        case {'Newton', 'regNewton'}

            x_exp_term_1 = x_1.*exp_term_1;
            alpha_x_exp_term_1 = alpha_1.*x_exp_term_1;
            sum_alpha_x_exp_term_1 = sum(alpha_x_exp_term_1);
            x_exp_term_2 = x_2.*exp_term_2;
            alpha_x_exp_term_2 = alpha_2.*x_exp_term_2;
            sum_alpha_x_exp_term_2 = sum(alpha_x_exp_term_2);            
           
            while abs(g_1) > tol_linesearch 
                
                nom_g_2_1 = sum(alpha_1.*alpha_x_exp_term_1) * sum_x_exp_term_1 - sum_alpha_x_exp_term_1^2;
                denom_g_2_1 = sum_x_exp_term_1^2;
                nom_g_2_2 = sum(alpha_2.*alpha_x_exp_term_2) * sum_x_exp_term_2 - sum_alpha_x_exp_term_2^2;
                denom_g_2_2 = sum_x_exp_term_2^2;            
                g_2 = nom_g_2_1 / denom_g_2_1 + nom_g_2_2 / denom_g_2_2; 
                
                if strcmp(solver.linesearch_optimizer, 'Newton')
                    t = t - g_1/g_2;
                else
                    t = t - g_1/(g_2 + solver.Newton_damping_factor * sqrt(abs(g_1)));
                end
    
                assert( abs(t) < solver.max_abs_t )
                
                if solver.stabilize_exp
                    exp_term_1 = exp(t*(alpha_1_max-alpha_1)); 
                else
                    exp_term_1 = exp(-t*alpha_1); 
                end
                x_exp_term_1 = x_1.*exp_term_1; 
                alpha_x_exp_term_1 = alpha_1.*x_exp_term_1; 
                sum_x_exp_term_1 = sum(x_exp_term_1); 
                sum_alpha_x_exp_term_1 = sum(alpha_x_exp_term_1); 
                if solver.stabilize_exp
                    exp_term_2 = exp(t*(alpha_2_max-alpha_2));
                else
                    exp_term_2 = exp(-t*alpha_2);
                end
                x_exp_term_2 = x_2.*exp_term_2; 
                alpha_x_exp_term_2 = alpha_2.*x_exp_term_2; 
                sum_x_exp_term_2 = sum(x_exp_term_2); 
                sum_alpha_x_exp_term_2 = sum(alpha_x_exp_term_2); 
                g_1 = beta - sum_alpha_x_exp_term_1 / sum_x_exp_term_1 - sum_alpha_x_exp_term_2 / sum_x_exp_term_2; 
                % debug
                %plot(t,eval_g(t),'*','color','red');
                %steps = steps + 1;
            end

        % debug
        %disp(steps)
        %fprintf('Stepsize: %f,   initial value (from eucl. proj.): %f \n', t, t_eucl)
        %plot(t,eval_g(t),'*','color','green');

        %if steps > 30
        %    fprintf('More than 30 steps of Newton`s method for linesearch (%d) \n', steps)
        %end
        
        t_opt = t;

     end
   

    
end