function J = costLanding2(xx)
    n  = (length(xx) - 2)/7;
    wf = xx(4*n);
    uf = xx(3*n);

    
    % J = wf^2;

    %
    wlim = 100*0.3048/60;
    ulim = 15*0.5144*0;

    e1 = wf - wlim;
    e2 = uf - ulim;
    J1  = e1^2 + e2^2;
    %}
    % J = 0;
    J2 = 10*(wf)^2 +uf^2;
    J=J2;
end