% Visualization of the step size found by product_probability_simplex_entropy_stepsize
function plot_product_probability_simplex_entropy_linesearch_function(x_1, x_2, alpha_1, alpha_2, beta)

    g = @(t) beta*t + log( sum( x_1.* exp(-t*alpha_1) ) )...
                                     + log( sum( x_2.* exp(-t*alpha_2) ) );
    
    t_SPS = ( alpha_1'*x_1 + alpha_2'*x_2 - beta)/ (norm(alpha_1)^2 + norm(alpha_2)^2);
    
    ts = linspace(t_SPS-10, t_SPS+10, 10);
    gs = zeros(size(ts));
    
    for i = 1:numel(ts)
        gs(i) = g(ts(i));
    end

    ts = ts(gs<Inf);
    gs = gs(gs<Inf);

    figure
    plot(ts, gs, 'color','blue')
    title('Linesearch objective')
    hold on

end