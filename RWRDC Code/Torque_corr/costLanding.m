function J = costLanding(xx)
    n  = 0.125*(length(xx) - 1);
    uf = xx(3*n);
    wf = xx(4*n);

    % Choose and modify whichever cost function seems suitable

    J1 = wf^2 + 0*uf^2;

    e1 = wf - 100*0.3048/60;
    e2 = uf - 0*0.5144;
    J2 = 100*e1^2 + e2^2;

    J = J2;
end