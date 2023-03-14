% uniform distribution on the probability simplex
function x = random_point_on_probability_simplex(dim)

    x = rand(dim, 1);
    x = -log(x);
    x = x / sum(x); 

end