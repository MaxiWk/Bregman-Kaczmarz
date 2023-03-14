% yields a random number i between 1 and length(p) with probability p(i)
function i = sampling(p)
    P = cumsum(p); 
    i = nnz(rand>P)+1;
end