clc; clear; close all;

% Initializing Parameters
rotor_init;

% Optimizer Parameters
opts = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12, ...
                    'Algorithm','trust-region-dogleg','Display','none', ...
                    'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);


% Finding Autorotation Trim of 6-DoF system to use as a guess
load trim_AR_guess_20.mat
gusart = out;
sts    = [0; 0; ufmin; v_i; omg_i]; % x_i, z_i don't matter

trARF = @(x, sts) state_eqns_6D(sts(1:4), [x; sts(5)], heli, par, 4);

[aTr,~,flg] = fsolve(@(x) trARF(x, sts), gusart, opts);

if any(flg == [1,2,3,4])
    fprintf('Autorotation Trim Solution Obtained\n')
end


% Load the solution from one of _ind files (Individual solutions).
[lgd, folderPath] = uigetfile('*.mat','Select mat files','Multiselect','off');
S = load(lgd);
% S = load('ThrPhs_MoI-1550_v-0_h-820_TAuend20kts.mat');

% Finding the number of phases
nd  = S.nd; 
nTD = S.nTD;
nP  = (length(S.t3D) - nTD)/nd;
D   = Dmat(nd - 1);

t3D  = S.t3D(nTD+1:end);
t3D  = t3D - t3D(1);
t3DM = reshape(t3D, nd, nP);
t3DM = t3DM - (ones(20,1) * t3DM(1,:));

% Main Loop
sol6DE = [];
t6DE   = [];
S1     = load('sol6DE1800_TRD_DC_poly.mat'); % For using as guess
% S1     = load('newgues.mat'); % For using as guess

for p = 1:nP
    fprintf('Phase-%d:\n', p)
    o=polyfit(t3DM(:,p),S.sol3D(nTD+(1:nd)+(p-1)*nd,end),6);
    k=polyder(o);
    l=polyder(k);
    z=polyval(l,t3DM(:,p));
    tV  = t3DM(:,p);
    sol = S.sol3D(nTD+(1:nd)+(p-1)*nd, :);
    D1  = 2*D/tV(end);

    attd      = zeros(nd, 8);   % Attitude matrix
    attd(:,1) = sol(:,8);       % thtP
    attd(:,4) = zeros(nd, 1);   % phiR
    attd(:,7) = zeros(nd, 1);   % psiY

    % Uncomment the following lines if a polynomial has to be fit on the
    % pitch angle solution.
    % ordPoly = 6; % Order of Polynomial
    % attd(:,1) = polyval(polyfit(tV, attd(:,1), 6), tV);

    attd(:,2) = polyval(k,t3DM(:,p));   % thtPD
    attd(:,3) =z;   % thtPD2
    attd(:,5) = D1   * attd(:,4);   % phiRD
    attd(:,6) = D1^2 * attd(:,4);   % phiRD2
    attd(:,8) = D1   * attd(:,7);   % psiYD

    % Extrapolating the end values to smoothen out any weird jumps
    attd([1,end],2) = interp1(tV(2:end-1), attd(2:end-1,2), tV([1,end]),'linear','extrap');
    attd([1,end],3) = interp1(tV(2:end-1), attd(2:end-1,3), tV([1,end]),'linear','extrap');
    attd([1,end],5) = interp1(tV(2:end-1), attd(2:end-1,5), tV([1,end]),'linear','extrap');
    attd([1,end],6) = interp1(tV(2:end-1), attd(2:end-1,6), tV([1,end]),'linear','extrap');

    solM = zeros(nd, 11);
    flgV = zeros(nd, 1);
    cnt  = 0;

    for i = 1:nd
        sts  = [sol(i,1:3), 0, sol(i,4:5), attd(i, [1,2,4,5])]';
        ctrl = [sol(i,7), attd(i, [3,6,8])]';

        if (p-1)*20+i ~= length(S1.sol6DE)
            gus = S1.sol6DE((p-1)*20+i ,[16:22, 11:14])';
        end

        while cnt < 3
            [out,~,flg] = fsolve(@(x) extend3D26D(x, sts, ctrl, heli, par), gus, opts);

            solM(i,:) = out;
            flgV(i)   = flg;

            cnt = cnt + 1;
            if any(flg == 1:4)
                gus = out;
                break
            elseif cnt == 2
                tmp = fsolve(@(x) trARF(x, [0;0; sol(i,3); 0; sol(i,5)]), gusart, opts);
                gus = tmp([2:8, 11:14]);
            elseif cnt == 3
                gus = aTr([2:8, 11:14]);
            end
        end

        fprintf('Node %2.0f/%2.0f Flag = %2.0f No. of tries = %d\n', i,nd, flg, cnt)
        cnt = 0;
    end

    fprintf('Phase-%d Completed!!! ',p)
    fail = ~ismember(flgV, 1:4);
    fprintf('No. of failed cases = %d\n\n', sum(fail))

    indF = setdiff(1:nd, find(fail));

    tmp = [sol(:,1:3), zeros(nd,1), sol(:,4:5), attd(:, [1,2,4,5]),...
           solM(:,8:11), sol(:,7), solM(:, 1:7)];
    sol6DE = [sol6DE; tmp(indF,:)];
    t6DE   = [t6DE; tV(indF) + S.t3D(nTD+1+(p-1)*nd)];
end


% If the solution has to be saved for plotting/other purposes, set the
% variable "toSave" as 1. Otherwise, the code will stop run without saving.
% Change filename in variable "fname" if required.
toSave = 1;

if toSave == 1
    % Name of the .mat file
    fname = 'sol6DE_1800.mat';

    save(fname, 'sol6DE', 't6DE', 'nd', 'nP')
else
    return
end
toSave2 = 0;

if toSave2 == 1
    % Name of the .mat file
    fname = 'sol6DE_g.mat';

    save(fname, 'sol6DE', 't6DE')
else
    return
end