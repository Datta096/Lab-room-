function [dx, op] = state_eqns_6D(sts, ctrl, heli, par, md)
    % Always States:
    % sts = [x, z, u, v]
    %
    % Always Controls:
    % ctrl = [tht0, thtc, thts, thtT, bet0, betc, bets, betT]
    %
    % Mode 1 -> Optimizer
    % sts  = 14x1 | [sts, w, omg, tht, thtdot, phi, phidot, 4xlmb]
    % ctrl =  8x1 | [ctrl]
    % dx   = 18x1 | [dxS, dxMb]
    % dxS  = 14x1 | [dot{x, z, u, v, w, omg, 4xlmb, tht, thtdot, phi, phidot}]
    % dxMb =  4x1 | [Mb0, Mbc, Mbs, MbT]
    %
    % Mode 2 -> Full Trim
    % sts  =  6x1 | [sts, w, omg]
    % ctrl = 14x1 | [ctrl, tht, phi, 4xlmb]
    % dx   = 14x1 | [dxS, dxM, dxMb]
    % dxS  =  7x1 | [dot{u, v, w, 4xlmb}]
    % dxM  =  3x1 | [Mx, My, Mz]
    % dxMb =  4x1 | [Mb0, Mbc, Mbs, MbT]
    %
    % Mode 3 -> Time March
    % sts  = 14x1 | [sts, w, omg, tht, thtdot, phi, phidot, 4xlmb]
    % ctrl =  8x1 | [ctrl]
    % dx   = 14x1 | [dxS]
    % dxS  = 14x1 | [dot{x, z, u, v, w, omg, 4xlmb, tht, thtdot, phi, phidot}]
    %
    % Mode 4 -> ROD Minimization
    % sts  =  4x1 | [sts]
    % ctrl = 16x1 | [ctrl, tht, phi, 4xlmb, w, omg]
    % dx   = 15x1 | [dxS, dxM, dxMb]
    % dxS  =  8x1 | [dot{u, v, w, omg, 4xlmb}]
    % dxM  =  3x1 | [Mx, My, Mz]
    % dxMb =  4x1 | [Mb0, Mbc, Mbs, MbT]
    %
    % Mode 5 -> Beta Equation
    % sts  = 14x1 | [sts, w, omg, tht, thtdot, phi, phidot, 4xlmb]
    % ctrl =  8x1 | [ctrl]
    % dx   =  4x1 | [dxMb]
    % dxMb =  4x1 | [Mb0, Mbc, Mbs, MbT]

    rotor = heli.rotor;
    hvTRp = heli.hvTRp;

    FWI = 0; % Flag for Wake Interaction
    FAW = 1; % Flag for Aero Wings

    % Always States
    u   = sts(3); % NED
    v   = sts(4);

    % Always Controls
    tht0 = ctrl(1);
    thtc = ctrl(2);
    thts = ctrl(3);
    thtT = ctrl(4);
    bet0 = ctrl(5) + rotor.precon; % Rotor pre-cone angle
    betc = ctrl(6);
    bets = ctrl(7);
    betT = ctrl(8);

    if any(md == [1,3,5])
        w      = sts(5);
        omg    = sts(6);
        ang    = sts(7);
        thtdot = sts(8);
        phiR   = sts(9);
        phidot = sts(10);
        lmb0   = sts(11);
        lmbs   = sts(12);
        lmbc   = sts(13);
        lmbT   = sts(14);
    elseif md == 2
        w    = sts(5);
        omg  = sts(6);
        ang  = ctrl(9);
        phiR = ctrl(10);
        lmb0 = ctrl(11);
        lmbs = ctrl(12);
        lmbc = ctrl(13);
        lmbT = ctrl(14);

        phidot = 0;
        thtdot = 0;
    elseif md == 4
        ang  = ctrl(9);
        phiR = ctrl(10);
        lmb0 = ctrl(11);
        lmbs = ctrl(12);
        lmbc = ctrl(13);
        lmbT = ctrl(14);
        w    = ctrl(15);
        omg  = ctrl(16);

        phidot = 0;
        thtdot = 0;
    end



    % Non-dimensional velocities w.r.t. the rotor
    angR = ang + rotor.alp_s;
    TpSp = omg*rotor.R;
    lmbd = (u*sin(angR) + w*cos(angR))/TpSp;
    mu   = (u*cos(angR) - w*sin(angR))/TpSp;

    % Orientation of Tip Path Plane
    alplong = angR - betc;
    alplat  = phiR + bets;
    den     = sqrt(betc^2 + bets^2 + 1);

    % Adding rotor velocity to u, w to be used for rest of the helicopter
    ur = lmb0*sin(alplong)*TpSp*FWI;
    wr = lmb0*cos(alplong)*TpSp*FWI;
    u  = u - ur;
    w  = w - wr;
    V  = sqrt(u^2 + v^2 + w^2);
    if V == 0
        beta = 0; % Side-slip
    else
        beta = asin(v/V);
    end
    % beta = 2*pi/180;
    gma  = atan2(w, u);
    alpF = ang + gma; % Fuselage Angle of attack

    % Computing Flap Harmonics
    kb  = rotor.kb;
    ehM = rotor.hngofs;
    nu  = kb*rotor.Nb/(rotor.MoI*omg^2) + 1.5*ehM/(1 - ehM);

    if heli.f1 == 0
        npsi    = 1;
        psi_arr = 0;
    else
        npsi    = rotor.npsi;
        psi_arr = linspace(0, 2*pi, npsi)';
    end

    r   = rotor.xk;
    wk  = rotor.wk;
    eMR = rotor.e;

    frcMat = zeros(npsi, 3);
    momMat = zeros(npsi, 3);
    othMat = zeros(npsi, 5);
    flpMat = zeros(npsi, 3);
    MdelMt = zeros(npsi, 3);

    for i = 1:npsi
        psi = psi_arr(i);

        % Flap and its Derivative
        bet    = bet0 + betc*cos(psi) + bets*sin(psi);
        betDot =      - betc*sin(psi) + bets*cos(psi);

        % Total inflow
        lmbi = lmb0 + lmbc*r*cos(psi) - lmbs*r*sin(psi);
        lmda = lmbi - lmbd;

        % Non-dimensionalised Velocities w.r.t. Blade Section
        UTbar = r + mu*sin(psi) - heli.zMR*(phidot*cos(psi) + thtdot*sin(psi))/TpSp;
        UPbar = lmda + r*betDot + mu*bet*cos(psi);
        UPbar = UPbar + (thtdot*(r*cos(psi) - heli.xMR/rotor.R) - phidot*r*sin(psi))/omg;
        U2bar = UPbar.^2 + UTbar.^2;

        % Blade pitch moment
        dM   = -0.07*(0.5*par.rho*U2bar*TpSp^2*rotor.c^2);
        Mdel = wk*dM * 0.5*(1 - eMR); % r-integrated

        MdelMt(i,1) = Mdel;
        MdelMt(i,2) = Mdel*cos(psi);
        MdelMt(i,3) = Mdel*sin(psi);
    end

    Idel = rotor.MB*rotor.c^2/48;
    Kdel = 30000; %25*omg^2*Idel;
    del0 = trapz(psi_arr, MdelMt(:,1))/Kdel;
    delc = trapz(psi_arr, MdelMt(:,2))/(pi*(Kdel - Idel*omg^2));
    dels = trapz(psi_arr, MdelMt(:,3))/(pi*(Kdel - Idel*omg^2));

    for i = 1:npsi
        psi = psi_arr(i);

        % Flap and its Derivative
        bet    = bet0 + betc*cos(psi) + bets*sin(psi);
        betDot =      - betc*sin(psi) + bets*cos(psi);

        % Total inflow
        lmbi = lmb0 + lmbc*r*cos(psi) - lmbs*r*sin(psi);
        lmda = lmbi - lmbd;

        % Non-dimensionalised Velocities w.r.t. Blade Section
        UTbar = r + mu*sin(psi) - heli.zMR*(phidot*cos(psi) + thtdot*sin(psi))/TpSp;
        UPbar = lmda + r*betDot + mu*bet*cos(psi);
        UPbar = UPbar + (thtdot*(r*cos(psi) - heli.xMR/rotor.R) - phidot*r*sin(psi))/omg;
        U2bar = UPbar.^2 + UTbar.^2;

        % Induced Angle of Attack
        phi = atan2(UPbar, UTbar);
        % Total Pitch Angle of the Blade Section
        tht = tht0 + thtc*cos(psi) + thts*sin(psi);
        % tht = tht + del0 + delc*cos(psi) + dels*sin(psi);
        % Angle of Attack of the Blade Section
        alfa = tht + rotor.tht_rt + rotor.twist*(r - eMR) - phi;

        if heli.f2 == 0
            Cl = rotor.Cla*alfa + rotor.Cl0;
        elseif heli.f2 == 1
            Cl = zeros(size(alfa));
            for j = 1:length(alfa)
                Cl(j) = getCl(alfa(j));
            end
        elseif heli.f2 == 2
            Cl = rotor.Cla*sin(alfa).*cos(alfa) + rotor.Cl0;
        end
        Cd = rotor.Cd0 + rotor.Cda*alfa + rotor.Cda2*alfa.^2;

        % Cl(Cl > 1.2) = 1.2;

        % [Cl, Cd] = getClCd(1, par, alfa*180*pi, 0.4);

        dCFz = Cl.*cos(phi) - Cd.*sin(phi);
        dCFx = Cl.*sin(phi) + Cd.*cos(phi);

        frcMat(i,1) = wk*(U2bar.* -dCFx) * 0.5*(1 - eMR)*sin(psi);
        frcMat(i,2) = wk*(U2bar.*  dCFx) * 0.5*(1 - eMR)*cos(psi);
        frcMat(i,3) = wk*(U2bar.*  dCFz) * 0.5*(1 - eMR);

        momMat(i,1) = wk*(U2bar.*  dCFz.*r) * 0.5*(1 - eMR)*sin(psi);
        momMat(i,2) = wk*(U2bar.* -dCFz.*r) * 0.5*(1 - eMR)*cos(psi);
        momMat(i,3) = wk*(U2bar.* -dCFx.*r) * 0.5*(1 - eMR);

        othMat(i,1) =  rotor.kb*bet*sin(psi)*rotor.Nb;
        othMat(i,2) = -rotor.kb*bet*cos(psi)*rotor.Nb;
        othMat(i,3) =  wk*(U2bar.^1.5 .* dCFx) * 0.5*(1 - eMR);

        % Blade pitch moment
        dM   = -0.07*(0.5*par.rho*U2bar*TpSp^2*rotor.c^2);
        Mdel = wk*dM * 0.5*(1 - eMR); % r-integrated
        othMat(i,2) = othMat(i,2) - Mdel*sin(psi);

        Mbeta       = 0.5*rotor.gma*wk*(U2bar.*  dCFz.*r) * 0.5*(1 - eMR)/rotor.Cla;
        flpMat(i,1) = bet0 + nu*bet - Mbeta + 0.5*rotor.MB*par.g*rotor.R/(rotor.MoI*omg^2);
    end

    frcMat      =  0.5*rotor.sigma * frcMat      * rotor.term*omg^2;
    momMat      =  0.5*rotor.sigma * momMat      * rotor.term*omg^2*rotor.R;
    othMat(:,3) =  0.5*rotor.sigma * othMat(:,3) * rotor.term*omg^3*rotor.R;
    othMat(:,4) =  frcMat(:,3)*ehM*rotor.R.*sin(psi_arr);
    othMat(:,5) = -frcMat(:,3)*ehM*rotor.R.*cos(psi_arr);
    flpMat(:,2) =  flpMat(:,1).*cos(psi_arr);
    flpMat(:,3) =  flpMat(:,1).*sin(psi_arr);

    if heli.f1 == 0
        HMR = frcMat(1);
        SMR = frcMat(2);
        TMR = frcMat(3);

        MxAeMR = momMat(1);
        MyAeMR = momMat(2);
        QMR    = momMat(3);

        MxS  = othMat(1);
        MyS  = othMat(2);
        PMR  = othMat(3);
        MxHO = othMat(4);
        MyHO = othMat(5);

        Mbet0 = flpMat(1);
        Mbetc = flpMat(2);
        Mbets = flpMat(3);
    else
        HMR = (0.5/pi)*trapz(psi_arr, frcMat(:,1));
        SMR = (0.5/pi)*trapz(psi_arr, frcMat(:,2));
        TMR = (0.5/pi)*trapz(psi_arr, frcMat(:,3));

        MxAeMR = (0.5/pi)*trapz(psi_arr, momMat(:,1));
        MyAeMR = (0.5/pi)*trapz(psi_arr, momMat(:,2));
        QMR    = (0.5/pi)*trapz(psi_arr, momMat(:,3));

        MxS  = (0.5/pi)*trapz(psi_arr, othMat(:,1));
        MyS  = (0.5/pi)*trapz(psi_arr, othMat(:,2));
        PMR  = (0.5/pi)*trapz(psi_arr, othMat(:,3));
        MxHO = (0.5/pi)*trapz(psi_arr, othMat(:,4));
        MyHO = (0.5/pi)*trapz(psi_arr, othMat(:,5));

        Mbet0 = (0.5/pi)*trapz(psi_arr, flpMat(:,1));
        Mbetc = (0.5/pi)*trapz(psi_arr, flpMat(:,2));
        Mbets = (0.5/pi)*trapz(psi_arr, flpMat(:,3));
    end

    rot = [cos(beta), -sin(beta); sin(beta), cos(beta)];
    tmp = rot*[HMR; SMR];
    HMR = tmp(1);
    SMR = tmp(2);

    tmp = rot*[MxS; MyS];
    MxS = tmp(1);
    MyS = tmp(2);

    tmp  = rot*[MxHO; MyHO];
    MxHO = tmp(1);
    MyHO = tmp(2);

    % Components of thrust along the coordinate axes
    TxMR =  TMR*(betc - rotor.alp_s)/den;
    TyMR =  TMR*bets/den;
    TzMR = -TMR/den;

    CTLM = [TMR*rotor.R; MxAeMR; -MyAeMR]/(rotor.term*omg^2*rotor.R);


    % Fuselage
    Lfus    = 0.5*par.rho*V^2*polyval(heli.CLfus, alpF*180/pi);
    Dfus    = 0.5*par.rho*V^2*polyval(heli.CDfus, alpF*180/pi);
    Pfus    = 0.5*par.rho*V^3*polyval(heli.CDfus, alpF*180/pi);


    % H-Tail
    if FAW == 1
        SH = hvTRp.SH;
    else
        SH = 0;
    end
    ClH   = hvTRp.ClaH*(ang + gma + hvTRp.alpH0);
    LH    = 0.5*par.rho*V^2*ClH*SH;
    DH    = 0.5*par.rho*V^2*(hvTRp.CD0H + hvTRp.KH*ClH^2)*SH;


    % V-Tail
    if FAW == 1
        SV = hvTRp.SV;
    else
        SV = 0;
    end
    ClV   = hvTRp.ClaV*(beta + hvTRp.alpV0);
    LV    = 0.5*par.rho*V^2*ClV*SV; % Assume to left for now
    DV    = 0.5*par.rho*V^2*(hvTRp.CD0V + hvTRp.KV*ClV^2)*SV;


    % Tail Rotor
    RTR   = hvTRp.RTR;
    omgT  = hvTRp.TRGR * omg;
    del3  = hvTRp.del3;
    eTR   = hvTRp.eTR;
    r     = hvTRp.xk;

    muTR   = (u*cos(ang) - w*sin(ang))/(omgT*RTR);
    frcMat = zeros(npsi, 3);
    othMat = zeros(npsi, 5);

    for i = 1:npsi
        psi = psi_arr(i);

        % Total inflow
        lmda = lmbT - v/(omgT*RTR);

        % Non-dimensionalised Velocities w.r.t. Blade Section
        UTbar = r + muTR*sin(psi);
        UPbar = lmda;
        U2bar = UPbar.^2 + UTbar.^2;

        % Induced Angle of Attack
        phi  = atan2(UPbar, UTbar);
        % Pitch angle
        tht  = thtT + hvTRp.thtRtTR + hvTRp.twsTR*(r - eTR);
        % Angle of attack
        alfa = tht - betT*del3 - phi;

        Cl = hvTRp.Cl0 + hvTRp.Cla*alfa;
        Cd = hvTRp.Cd0 + hvTRp.Cda*alfa + hvTRp.Cda2*alfa.^2;

        dCFz = Cl.*cos(phi) - Cd.*sin(phi);
        dCFx = Cl.*sin(phi) + Cd.*cos(phi);

        frcMat(i,1) = wk*(U2bar.* -dCFx) * 0.5*(1 - eTR)*sin(psi);
        frcMat(i,2) = wk*(U2bar.*  dCFx) * 0.5*(1 - eTR)*cos(psi);
        frcMat(i,3) = wk*(U2bar.*  dCFz) * 0.5*(1 - eTR);

        othMat(i,1) = wk*(U2bar.*       dCFz.*r) * 0.5*(1 - eTR)*sin(psi);
        othMat(i,2) = wk*(U2bar.*      -dCFz.*r) * 0.5*(1 - eTR)*cos(psi);
        othMat(i,3) = wk*(U2bar.*      -dCFx.*r) * 0.5*(1 - eTR);
        othMat(i,4) = wk*(U2bar.^1.5 .* dCFx)    * 0.5*(1 - eTR);

        Mbeta       = 0.5*hvTRp.gmaTR*wk*(U2bar.*  dCFz.*r) * 0.5*(1 - eTR)/hvTRp.Cla;
        othMat(i,5) = betT - Mbeta;
    end

    frcMat        = 0.5*hvTRp.sigTR * frcMat        * (par.rho*pi*RTR^4)*omgT^2;
    othMat(:,1:3) = 0.5*hvTRp.sigTR * othMat(:,1:3) * (par.rho*pi*RTR^5)*omgT^2;
    othMat(:,4)   = 0.5*hvTRp.sigTR * othMat(:,4)   * (par.rho*pi*RTR^5)*omgT^3;

    if heli.f1 == 0
        HTR = frcMat(1);
        STR = frcMat(2);
        TTR = frcMat(3);

        MxTR =  othMat(1);
        MyTR =  othMat(3);
        MzTR = -othMat(2);
        PTR  =  othMat(4);
        MbtT =  othMat(5);
    else
        HTR = (0.5/pi)*trapz(psi_arr, frcMat(:,1));
        STR = (0.5/pi)*trapz(psi_arr, frcMat(:,2));
        TTR = (0.5/pi)*trapz(psi_arr, frcMat(:,3));

        MxTR = (0.5/pi)*trapz(psi_arr,  othMat(:,1));
        MyTR = (0.5/pi)*trapz(psi_arr,  othMat(:,3));
        MzTR = (0.5/pi)*trapz(psi_arr, -othMat(:,2));
        PTR  = (0.5/pi)*trapz(psi_arr,  othMat(:,4));
        MbtT = (0.5/pi)*trapz(psi_arr,  othMat(:,5));
    end

    CTTR = TTR/(par.rho*pi*RTR^4*omgT^2);
    VTTR = sqrt(lmbT^2 + muTR^2);


    % Dynamic Inflow Model
    lmda    = lmb0 - lmbd;
    lmbh    = sqrt(0.5*abs(CTLM(1)));
    lmb_bar = lmda/lmbh;

    % if mu<(0.1*0.5144/TpSp)
    %     mu_bar = (0.1*0.5144/TpSp)/lmbh;
    % else
        mu_bar  = mu/lmbh;
    % end

    %{
    if (lmb_bar >= -1) && (lmb_bar <= 0.6378)
        g_lmda = 1/(2 + lmb_bar)^2 - lmb_bar^2 + ...
                 (1 + lmb_bar)* (0.109 + 0.217*(lmb_bar - 0.15)^2);
        g_prime = (-2/(2+lmb_bar)^3)+0.049-1.696*lmb_bar+0.651*lmb_bar^2;
    else
        g_lmda = 0;
        g_prime = 0;
    end

    if (mu_bar >= 0) && (mu_bar <= 0.707)
    % if (mu_bar^2 >= 0) && (mu_bar^2 <= 0.5)
        f_mu = 1 - 2*mu_bar^2;
        f_prime = -4*mu_bar;
        % f_mu = 0;
    elseif mu_bar < 0
        f_mu = 1;
        f_prime = 0;
    else
        f_mu = 0;
        f_prime = 0;
    end
    %}

    %
    g_lmda  = 0.40223*exp(-3.081*(lmb_bar + 0.15016).^2);
    g_prime = -2.47853*(lmb_bar + 0.15016).*exp(-3.081*(lmb_bar + 0.15016).^2);
    
    if mu_bar < 0
        f_mu    = 1;
        f_prime = 0;
    else
        f_mu    = exp(-3.45768*mu_bar.^2);
        f_prime = -6.91536*mu_bar.*exp(-3.45768*mu_bar.^2);
    end
    %}


    % if (mu_bar >= 0) && (mu_bar <= 1/sqrt(3))
    %     f_mu = 1 - 3*mu_bar^2;
    % else
    %     f_mu = 0;
    % end
    % f_mu = 0; g_lmda = 0;
    
    VT   = sqrt(lmda^2 + mu^2 + lmbh^2*g_lmda*f_mu);
    % VM   = (mu^2 + lmda*(lmda + lmb0))/VT;

    VM = (mu^2 + lmda*(lmda + lmb0) + lmbh^2*f_mu*(g_lmda+0.5*g_prime*(lmb0/lmbh)))/...
         (VT - 0.5*lmb0*g_lmda*f_mu + 0.25*lmb0*(lmb_bar*g_prime*f_mu + mu_bar*g_lmda*f_prime));

    Vmat = diag([VT, VM, VM]);
    % sx   = sin(atan2(lmda, mu));
    sx = abs(lmda/norm([lmda, mu]));
    % sx = lmda/norm([lmda, mu]);
    % sx = 1;
    sxp  = 1/(1 + sx);
    sxm  = 1 - sx;
    Mmat = heli.M;
    t    = heli.t;
    Lmat = [                 0.5,     0, 0.25*t*sqrt(sxm*sxp);
                               0, -4*sxp,                     0;
            0.25*t*sqrt(sxm*sxp),     0,              4*sx*sxp];


    % Output Structure
    op.MR.frc = [HMR, SMR, TMR, HMR + TxMR, SMR + TyMR, TzMR];
    op.MR.mom = [MxAeMR, MyAeMR, QMR, MxS, MyS, MxHO, MyHO];
    op.MR.oth = [PMR, bet0, betc, bets];

    op.TR.frc = [HTR, STR, TTR];
    op.TR.oth = [PTR, MxTR, MyTR, MzTR];

    op.FS.all = [Lfus, Dfus, Pfus];

    op.CS.htl = [LH, DH];
    op.CS.vtl = [LV, DV];

    op.oth = [PMR + PTR + Pfus, alpF*180/pi, lmb_bar, mu_bar, lmbh, lmbd, lmda,...
              del0, delc, dels];

    % Derivatives
    zdot = w;

    xdot = u;
    
    udot = -TMR*sin(alplong) + HMR*cos(alplong) + HTR*cos(ang) - STR*sin(ang);
    udot = (udot + (Lfus + LH)*sin(gma) - (Dfus + DH + DV)*cos(gma))/heli.MF;

    vdot = (TMR*sin(alplat) + SMR*cos(alplat) - LV*cos(phiR) - TTR*cos(phiR))/heli.MF;

    wdot = TMR*cos(alplong) + HMR*sin(alplong) + HTR*sin(ang) + STR*cos(ang);
    wdot = par.g - (wdot + (Lfus + LH)*cos(gma) + (Dfus + DH + DV)*sin(gma))/heli.MF;

    CTLM(2:3) = -CTLM(2:3);
    lmbidot = omg*(Mmat\(CTLM - Vmat*(Lmat\[lmb0; lmbs; lmbc])));
    % lmbidot = omg*(Mmat\(CTLM - Vmat*(Lmat\[lmb0; 0; 0])));
    lmbTdot = omgT*(CTTR - 2*VTTR*lmbT)/Mmat(1,1);

    odot = QMR/rotor.MoI;

    % Moment Equations
    dxM = zeros(3,1);

    % Main Rotor
    MxFMR = (SMR + TyMR)*heli.zMR;
    MyFMR = -TzMR*heli.xMR - (HMR + TxMR)*heli.zMR;
    MzMR  = QMR + SMR*heli.xMR;

    % Fuselage
    Myfus = 0.5*par.rho*V^2*polyval(heli.CMfus, alpF*180/pi);
    Myfus = Myfus + (Lfus*cos(alpF) + Dfus*sin(alpF))*heli.xFS;

    % H-Tail
    MyH = (DH*cos(alpF) - LH*sin(alpF))*heli.zHT - (LH*cos(alpF) + DH*sin(alpF))*heli.xHT;

    % Tail Rotor
    MxTR2 = MxTR - TTR*heli.zTR;
    MyTR2 = MyTR - HTR*heli.zTR - STR*heli.xTR;
    MzTR2 = MzTR + TTR*heli.xTR + HTR*heli.yTR;

    % V-Tail
    MyV = DV*(cos(alpF)*heli.zVT - sin(alpF)*heli.xVT);
    MzV = LV*heli.xVT;

    % Lateral
    MxLat = -TzMR*heli.yMR - (Lfus*cos(alpF) + Dfus*sin(alpF))*heli.yFS - ...
            (LH*cos(alpF) + DH*sin(alpF))*heli.yHT - DV*sin(alpF)*heli.yVT;
    MzLat = -(SMR + TxMR)*heli.yMR + (Dfus*cos(alpF) - Lfus*sin(alpF))*heli.yFS +...
            (DH*cos(alpF) - LH*sin(alpF))*heli.yHT + DV*cos(alpF)*heli.yVT;

    dxM(1) = MxS + MxHO + MxFMR + MxTR2 + MxLat;
    dxM(2) = Myfus + MyS + MyHO + MyFMR + MyH + MyTR2 + MyV;
    dxM(3) = MzMR + MzTR2 + MzV + MzLat;

    % Adding to output
    op.MR.mom(8:10) = [MxFMR, MyFMR, MzMR];
    op.TR.oth(5: 7) = [MxTR2, MyTR2, MzTR2];
    op.FS.all(4)    = Myfus;
    op.CS.htl(3)    = MyH;
    op.CS.vtl(3)    = MyV;
    op.CS.vtl(4)    = MzV;
    op.oth(11)      = dxM(3);

    % State Equations
    if any(md == [1,3])
        thtddot = dxM(2)/heli.Iy;
        phiddot = dxM(1)/heli.Ix;

        % dxS = [xdot; zdot; udot; vdot; wdot; 0; lmbidot(1);0;0; lmbTdot; thtdot; thtddot; phidot; phiddot];
        dxS = [xdot; zdot; udot; vdot; wdot; odot; lmbidot; lmbTdot; thtdot; thtddot; phidot; phiddot];
        % dxS = [xdot; zdot; udot; vdot; wdot; 0; lmbidot; lmbTdot; thtdot; thtddot; phidot; phiddot];
    elseif md == 2
        dxS = [udot; vdot; wdot; lmbidot; lmbTdot];
    elseif md == 4
        dxS = [udot; vdot; wdot; odot; lmbidot; lmbTdot];
    else
        dxS = [];
    end

    % Moment Equations
    if any(md == [1,3,5])
        dxM = [];
    end
    
    % Flap Equations
    if any(md == [1,2,4,5])
        dxMb = [Mbet0; Mbetc; Mbets; MbtT];
    else
        dxMb = [];
    end
 
    % Creating the final output vector
    dx = [dxS; dxM; dxMb];

    % keyboard
end