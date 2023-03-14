% Visualization of the step size found by probability_simplex_entropy_stepsize.m
function plot_probability_simplex_entropy_linesearch_function(x, alpha, beta)

    x_alpha = x.*alpha;

    g = @(t) beta*t + log( sum( x.* exp(-t*alpha) ) );
    eval_g_1 = @(t) beta - sum(x_alpha.*exp(-t*alpha)) / sum(x.*exp(-t*alpha));
    
    t_SPS = ( alpha'*x - beta)/ norm(alpha)^2; 
    
    ts = linspace(t_SPS-2, t_SPS+2, 100);
    gs = zeros(size(ts));
    g1s = zeros(size(ts));
    
    for i = 1:numel(ts)
        gs(i) = g(ts(i));
    end
    
    for i = 1:numel(ts)
        g1s(i) = eval_g_1(ts(i));
    end
   
    
    plot(ts, g1s, 'color','blue')
    title('Derivative of linesearch objective')

    figure
    plot(ts, gs, 'color','blue')
    title('Linesearch objective')
    hold on
    
    plot(t_SPS, g(t_SPS), '*', 'color', 'red')

end