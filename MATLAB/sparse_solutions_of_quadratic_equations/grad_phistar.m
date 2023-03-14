% soft shrinkage
function x = grad_phistar(xstar, lambda)  
    x = max(abs(xstar)-lambda, 0).*sign(xstar); 
end


