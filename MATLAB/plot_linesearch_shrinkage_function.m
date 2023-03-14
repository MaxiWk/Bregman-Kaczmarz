% Visualization of the step size found by linesearch_shrinkage.m
function plot_linesearch_shrinkage_function(x,x_star,a,y,c)

    [~,~,t_opt] = linesearch_shrinkage(x,x_star,a,y,c);
    
    t_min = t_opt - 2;
    t_max = t_opt + 2;
    
    ts = linspace(t_min,t_max,100);
    
    obj = @(t) 0.5*norm(soft_shrinkage(c,x_star-a*t))^2 + y*t;

    obj_vals = zeros(size(ts));


    for i = 1:length(ts)
        obj_vals(i) = obj(ts(i));
    end

    figure
    
    plot(ts,obj_vals)
    hold on
    plot(t_opt, obj(t_opt),'*','color','red');

    function v = soft_shrinkage(param, u)
        v = sign(u) .* max(abs(u)-param, 0);
    end
    
end