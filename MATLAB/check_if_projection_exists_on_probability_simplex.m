% checks if the condition in line 8 of Algorithm 1 is fulfilled in the
% case of the probability simplex constraint
function bool = check_if_projection_exists_on_probability_simplex(alpha, beta)
    tol = 1e-15;
    if all(abs(alpha-beta) < tol)
        bool = true;
        return;
    elseif any(alpha-beta > 0) && any(alpha-beta < 0)
        bool = true;
        return
    else
        bool = false;
    end
end