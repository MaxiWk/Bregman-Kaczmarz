% right inverse of soft shrinkage (grad_phistar) with least absolute value
function xstar = grad_phistar_inv(x, lambda)
    xstar = zeros(size(x));
    xstar(x>eps) = x(x>eps) + lambda;
    xstar(x<-eps) = x(x<eps) - lambda;
end