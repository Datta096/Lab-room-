function dx = LDM_state_eqns(sts, ctrl, heli, par, md)
    % sts  = [x, z, u, w, omg]
    % ctrl = [tht0, thtP]
    %
    % Mode 1 -> Optimizer / Time March
    % sts  = 6x1 | [sts, lmbi]
    % ctrl = 2x1 | [ctrl]
    % dx   = [dot{x, z, u, w, omg, lmbi}]
    %
    % Mode 2 -> Trim
    % sts  = 5x1 | [sts]
    % ctrl = 3x1 | [ctrl, lmbi]
    % dx   = [dot{u, w, lmbi}]
    %
    % Mode 3 -> Autorotation Trim
    % sts  = 4x1 | [x, z, u, omg]
    % ctrl = 4x1 | [ctrl, lmbi, w]
    % dx   = [dot{u, w, lmbi, omg}]

    rotor = heli.rotor;

    u   = sts(3);
    if any(md == 1:2)
        w   = sts(4);
        omg = sts(5);
    elseif md == 3
        w   = ctrl(4);
        omg = sts(4);
    end

    tht0 = ctrl(1);
    thtP = ctrl(2);
    thtc = par.thtc;
    thts = par.thts;
    angR = thtP + rotor.alp_s;

    if md == 1
        lmbi = sts(6);
    elseif any(md == 2:3)
        lmbi = ctrl(3);
    end

    TpSp = omg*rotor.R;
    lmbd = (u*sin(angR) + w*cos(angR))/TpSp;
    mu   = (u*cos(angR) - w*sin(angR))/TpSp;

    lmda = lmbi - lmbd;

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

    resM = zeros(npsi, 3);

    for i = 1:npsi
        psi = psi_arr(i);

        UTbar = r + mu*sin(psi);
        U2bar = lmda^2 + UTbar.^2; % UPbar^2 + UTbar^2
	    phi   = atan2(lmda, UTbar);
        tht   = tht0 + thtc*cos(psi) + thts*sin(psi);
	    alfa  = tht + rotor.tht_rt + rotor.twist*(r - eMR) - phi;

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
        Cl(alfa >= 11*pi/180) = interp1([11,45]*pi/180, [rotor.Cla*11* pi/180, 0]+rotor.Cl0, alfa(alfa >= 11*pi/180));

        Cd = rotor.Cd0 + rotor.Cda*alfa + rotor.Cda2*alfa.^2;

        % Cl(Cl > 1.2) = 1.2;

        dCFz = Cl.*cos(phi) - Cd.*sin(phi);
        dCFx = Cl.*sin(phi) + Cd.*cos(phi);
        
    	resM(i,1) = wk*(U2bar.*  dCFz)    * 0.5*(1 - eMR);
    	resM(i,2) = wk*(U2bar.* -dCFx.*r) * 0.5*(1 - eMR);
        resM(i,3) = wk*(U2bar.* -dCFx)    * 0.5*(1 - eMR)*sin(psi);
    end

    resM      = 0.5*rotor.sigma * resM * rotor.term*omg^2;
    resM(:,2) = resM(:,2) * rotor.R;

    if heli.f1 == 0
        T = resM(1);
        Q = resM(2);
        H = resM(3);
    else
        T = (0.5/pi)*trapz(psi_arr, resM(:,1));
        Q = (0.5/pi)*trapz(psi_arr, resM(:,2));
        
        H = (0.5/pi)*trapz(psi_arr, resM(:,3));
    end
	
    CT      = T/(rotor.term*omg^2);
    lmbh    = sqrt(0.5*abs(CT));
    lmb_bar = lmda/lmbh;
    mu_bar  = mu/lmbh;

    if heli.f3 == 1
        % New Fit on the Original Peters-He Model
        g_lmda  = 0.40223*exp(-3.081*(lmb_bar + 0.15016).^2);
        % g_prime = -2.47853*(lmb_bar + 0.15016).*exp(-3.081*(lmb_bar + 0.15016).^2);
        
        if mu_bar < 0
            f_mu    = 1;
            % f_prime = 0;
        else
            f_mu    = exp(-3.45768*mu_bar.^2);
            % f_prime = -6.91536*mu_bar.*exp(-3.45768*mu_bar.^2);
        end
    else
        % Original Peters-He Model
        if (lmb_bar >= -1) && (lmb_bar <= 0.6378)
            g_lmda = 1/(2 + lmb_bar)^2 - lmb_bar^2 + ...
                (1 + lmb_bar)* (0.109 + 0.217*(lmb_bar - 0.15)^2);
        else
            g_lmda = 0;
        end

        if (mu_bar >= 0) && (mu_bar <= 0.707)
            f_mu = 1 - 2*mu_bar^2;
        elseif mu_bar < 0
            f_mu = 1;
        else
            f_mu = 0;
        end
    end

    VT = sqrt(lmda^2 + mu^2 + lmbh^2*g_lmda*f_mu);

    gma  = atan2(w,u);
    alpF = thtP + gma;
    V    = sqrt(u^2 + w^2);
    Dfus = 0.5*par.rho*V^2*polyval(heli.CDfus, alpF*180/pi);
    Lfus = 0.5*par.rho*V^2*(polyval(heli.CLfus, alpF*180/pi) + 1.152);
    Ma   = 128/(75*pi);

    udot    = -(T*sin(angR) + Dfus*cos(gma) - H*cos(angR) - Lfus*sin(gma))/heli.MF;
    wdot    = par.g - (Dfus*sin(gma) + T*cos(angR) - H*sin(angR) + Lfus*cos(gma))/heli.MF;
    lmbidot = omg*(CT - 2*VT*lmbi)/Ma;
    %  if omg > 38.43
    %     corr = polyval(par.p, 38.43);
    % elseif omg < 0.65*38.43
    %     corr = polyval(par.p, 0.65*38.43);
    % else
    %     corr = polyval(par.p, omg);
    % end
    omgdot  = (Q *0.9)/rotor.MoI;

  

    if md == 1
        xdot   = u;
        zdot   = w;

        dx = [xdot; zdot; udot; wdot; omgdot; lmbidot];
    elseif md == 2
        dx = [udot; wdot; lmbidot];
    elseif md == 3
        dx = [udot; wdot; lmbidot; omgdot];
    end
end