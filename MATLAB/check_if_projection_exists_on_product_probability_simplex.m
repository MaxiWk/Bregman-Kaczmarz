% checks if the condition in line 8 of Algorithm 1 is fulfilled in the case
% of the product of the probability simplex constraint with itself
function bool = check_if_projection_exists_on_product_probability_simplex(alpha_1, alpha_2, beta)

    tol = 1e-15;

    alphas = {alpha_1, alpha_2};

    % none of the alphas is constant and there is common point in the 
    % intervals (beta - {min,max}(alpha_i)) and ({min,max}(alpha(j))
    % for (i,j) in {(1,2), (2,1)} 
    for order = [[1,2]; [2,1]]
        alpha_i = alphas{order(1)};
        alpha_j = alphas{order(2)};
        if check_if_intervals_intersect([beta-max(alpha_i), beta-min(alpha_i)],...
                                            [min(alpha_j), max(alpha_j)])
            bool = true; return
        end
    end
    bool = false;

    % exactly one of the alphas is constant (alpha_i = c) and 
    % min(alpha_j) < beta - c < max(alpha_j) 
    for order = [[1,2]; [2,1]]
        alpha_i = alphas{order(1)};
        alpha_j = alphas{order(2)};
        if all(abs(alpha_i-alpha_i(1)) < tol) 
           d = beta - alpha_i(1);
           bool = (min(alpha_j) < d && d < max(alpha_j)); return
        end
    end

    % both alphas are constant (c_1,c_2) and c_1 + c_2 = beta
    if all( abs(alpha_1 - alpha_1(1)) < tol & abs(alpha_2 - alpha_2(1)) < tol)
        bool = (abs(alpha_1(1) + alpha_2(1) - beta) < tol); return 
    end

end


