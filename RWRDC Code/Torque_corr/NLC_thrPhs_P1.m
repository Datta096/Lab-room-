function [C, Ceq] = NLC_thrPhs_P1(xx, D, n, BC, heli, par)
    x    = xx(    1:  n);
    z    = xx(  n+1:2*n);
    u    = xx(2*n+1:3*n);
    w    = xx(3*n+1:4*n);
    omg  = xx(4*n+1:5*n);
    lmbi = xx(5*n+1:6*n);
    thtP = xx(6*n+1:7*n);
    % tf   = xx(end-2);
    tf   = xx(end);
    
    % tf here is from start of optimizer, i.e., after delay and collective
    % dump
    D1   = 2*D/tf;
    Xdot = zeros(6, n);
    sts  = reshape(xx(1:6*n), n, 6);

    tht0 = xx(end-1)*ones(n,1);
    tsp  = (1 +chgslb(n))*0.5*tf;
    % tht0(tsp<xx(end)) = xx(end-1);
    ctrl = [tht0, thtP];

    for i = 1:n
        Xdot(:,i) = LDM_state_eqnsta(sts(i,:), ctrl(i,:), heli, par, 1,(tsp(i)+5.5));
    end
    C            = zeros(4*n+3, 1);
    C(  n+1:2*n) = abs(D1*thtP) - heli.rotor.alpdotmax;
    C(    1:  n) = abs(D1*tht0) - heli.rotor.thtdotmax;
    C(2*n+1:3*n) =  (D1*D1*thtP) - 0.45;
    C(3*n+1:4*n) = -(D1*D1*thtP) - 0.6;
    C(4*n+1)       =  xx(end) - tf;
    % C(4*n+2)     = -u(n)+BC(9);
    % C(end-1)    = -u(n)+30;
    % C(end)      = w(n)-30;
    Ceq            = zeros(6*n, 1);
    Ceq(    1:  n) = D1*x    - Xdot(1, :)';
    Ceq(  n+1:2*n) = D1*z    - Xdot(2, :)';
    Ceq(2*n+1:3*n) = D1*u    - Xdot(3, :)';
    Ceq(3*n+1:4*n) = D1*w    - Xdot(4, :)';
    Ceq(4*n+1:5*n) = D1*omg  - Xdot(5, :)';
    Ceq(5*n+1:6*n) = D1*lmbi - Xdot(6, :)';
    Ceq(1:n:5*n+1) = 0;
    % Initial conditions
    Ceq2(1) = x(1)    - BC(1);
    Ceq2(2) = z(1)    - BC(2);
    Ceq2(3) = u(1)    - BC(3);
    Ceq2(4) = w(1)    - BC(4);
    Ceq2(5) = omg(1)  - BC(5);
    Ceq2(6) = lmbi(1) - BC(6);
    Ceq2(8) = thtP(1) - BC(8);

    % % Terminal conditions
    Ceq2( 9) = u(n)    - BC( 9);
    % Ceq2(10) = w(n)    - BC(10);
    Ceq2(11) = omg(n)  - BC(11);
    % Ceq2(12) = lmbi(n) - BC(12);
    % Ceq2(13) = tht0(n) - BC(13);
    % Ceq2(14) = thtP(n) - BC(14);
    Ceq = [Ceq; Ceq2(:)];
end