function J = costSSAR(xx)
    hi = xx(21);
    hf = xx(40);
    J  = (-hi + hf);
end