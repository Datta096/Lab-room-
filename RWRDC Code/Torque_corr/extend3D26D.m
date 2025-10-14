function out = extend3D26D(x, sts, ctrl, heli, par)
    % sts  = [x, z, u, v, w, omg, 2xthtP, 2xphiR]
    % x    = [thtc, thts, thtT, 4xbet, 4xlmb]
    % ctrl = [tht0, thtPD2, phiRD2, psiYD2]
    % out  = [eqn{Mx, My, Mz}, 4xMb, dot{4xlmb}]

    tPD2 = ctrl(2);
    pRD2 = ctrl(3);
    pYD2 = ctrl(4);
    ctrl = [ctrl(1); x(1:7)];
    sts  = [sts; x(8:11)];

    % sts  = [x, z, u, v, w, omg, 2xthtP, 2xphiR, 4xlmb]
    % ctrl = [4xtht, 4xbet]
    % out  = [dot{x, z, u, v, w, omg, 4xlmb, 2xthtP, 2xphiR}, 4xMb]
    [out, op] = state_eqns_6D(sts, ctrl, heli, par, 1);
    out = [out([14,12]); op.oth(11); out([15:18, 7:10])];
    out(1) = out(1) - pRD2;
    out(2) = out(2) - tPD2;
    out(3) = out(3)/heli.Iz - pYD2; % Scaled to attitude ddots
end