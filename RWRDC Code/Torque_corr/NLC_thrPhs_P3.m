function [C, Ceq] = NLC_thrPhs_P3(xx, D, n, BC, heli, par)
     x    = xx(    1:  n);
    z    = xx(  n+1:2*n);
    u    = xx(2*n+1:3*n);
    w    = xx(3*n+1:4*n);
    omg  = xx(4*n+1:5*n);
    lmbi = xx(5*n+1:6*n);
    tht0 = xx(6*n+1:7*n);
    thtP = xx(7*n+1:8*n);
    tf   = xx(end);
    
    % tf here is from start of optimizer, i.e., after delay and collective
    % dump.
    D1   = 2*D/tf;
    Xdot = zeros(6, n);
    sts  = reshape(xx(    1:6*n), n, 6);
    ctrl = reshape(xx(6*n+1:8*n), n, 2);
   

    for i = 1:n
        Xdot(:,i) = LDM_state_eqns(sts(i,:), ctrl(i,:), heli, par, 1);
    end

    C          = zeros(2*n, 1);
    C(n+1:2*n) = abs(D1*thtP) - heli.rotor.alpdotmax;
    C(  1:  n) = abs(D1*tht0) - heli.rotor.thtdotmax;
    C          = [C; 1 + z + (163.36*thtP*180/pi + 170.89)*1e-3]; % For Tail-Guard
    C          = [C; z + (74.956*thtP*180/pi + 1502.6)*1e-3]; % For Aft-Skid

    Ceq            = zeros(6*n, 1);
    Ceq(    1:  n) = D1*x    - Xdot(1, :)';
    Ceq(  n+1:2*n) = D1*z    - Xdot(2, :)';
    Ceq(2*n+1:3*n) = D1*u    - Xdot(3, :)';
    Ceq(3*n+1:4*n) = D1*w    - Xdot(4, :)';
    Ceq(4*n+1:5*n) = D1*omg  - Xdot(5, :)';
    Ceq(5*n+1:6*n) = D1*lmbi - Xdot(6, :)';

    % Initial conditions
    Ceq2 = xx((0:7)*n+1) - BC(1:end-1);
    % Terminal conditions
    % Ceq3 = z(n) - BC(end);
    Ceq3 = z(n) + (74.956*thtP(n)*180/pi + 1502.6)*1e-3 + BC(end); % For Aft-Skid
    Ceq  = [Ceq; Ceq2; Ceq3];
end