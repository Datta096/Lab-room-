  clc; clear; close all;

% Initializing Parameters
rotor_init;

% Optimizer Parameters
opts{1} = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',2e7,...
                       'MaxIterations',1e5,'Algorithm','sqp','StepTolerance',1e-8, ...
                       'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-6); 

opts{2} = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12, ...
                       'Algorithm','trust-region-dogleg','Display','none', ...
                       'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);

opts{3} = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',2e7,...
                       'MaxIterations',1e5,'Algorithm','interior-point','StepTolerance',1e-7                                                                                                                                                  , ...
                       'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-6);

odeopt  = odeset('RelTol',1e-10, 'AbsTol',1e-10);

z_i =1000; % Set the initial altitude here, in feet
u_i = 0;   % Set the initial forward velocity here, in kts
% Rest of the states will be taken automatically from rotor_init

%% Read Initializations for Optimizer

load delM_fit_1550.mat 
par.p = -9.3127e2; 
% par.p =0;
par.TD = [2.5,2];

% Finding Trim Solutions
sts    = [x_i; -z_i/m2ft; u_i*kts2ms; w_i; omg_i];
gusfsl = [0.1; 0.1; 0.001];

[hTr,~,flg] = fsolve(@(ctrl) LDM_state_eqns(sts, ctrl, heli, par, 2), gusfsl, opts{2});

if any(flg == 1:4)
    fprintf('Trim Solution Obtained...\n')
end
X0     = [sts; hTr(3)];
sts    = [x_i; -z_i/m2ft; ufmin; omg_i];
gusfsl = [gusfsl; 0.1];

[aTr,~,flg] = fsolve(@(ctrl) LDM_state_eqns(sts, ctrl, heli, par, 3), gusfsl, opts{2});
if any(flg == 1:4)
    fprintf('Autorotation Trim Solution Obtained...\n')
end


% Time Marching for Pilot Delay & Collective Dump. Set in rotor_init.
if all(par.TD == 0)
    % No delay
    BC1 = [X0; hTr(1:2); ufmin; aTr(4); omg_i; aTr([3,1,2])];
else
    if par.TD(1) ~= 0
        tspan = linspace(0, par.TD(1), 10);
        [t1, Xs1] = ode45(@(t,x) LDM_state_eqnsta(x, hTr(1:2), heli, par, 1,t), tspan, X0, odeopt);
        X0 = Xs1(end,:)';
    else
        t1 = []; Xs1 = [];
    end
    if par.TD(2) ~= 0
        tspan = linspace(0, par.TD(2), 10);
        [t2, Xs2] = ode45(@(t,x) LDM_state_eqnsta(x, [-8.4*d2r; hTr(2)], heli, par, 1,t+par.TD(1)), tspan, X0, odeopt);
    else
        t2 = []; Xs2 = [];
    end
    tspan = [t1(:); t2(:) + par.TD(1)];
    Xsol  = [Xs1, ones(length(t1),1)*[  hTr(1), hTr(2)];
             Xs2, ones(length(t2),1)*[-8.4*d2r, hTr(2)]];
    BC1   = [Xsol(end,:)'; ufmin; aTr(4); omg_i; aTr([3,1,2])];
end


% load newGsol3P.mat
% load GsollinTa.mat
load Gsoltd3.5_1.6.mat
% load Guess3P100ss2.mat
% load new3PG3.mat
% load new3PG4.mat % 103%, 2500ft solution

fprintf('Phase-1... ')
lb(4*nd+1:5*nd) = 0.85*omg_i*ones(nd,1);
lb1 = lb; ub1 = ub;
lb1(4*nd+1:5*nd) = min(BC1(4), 0.85*omg_i);
Gsol1P= [Gsol1P(1:end-1);0*d2r;Gsol1P(end)];
% Gsol1P(end)   = 7;
ub1(4*nd+1:5*nd) = 1.01*omg_i;
lb1(6*nd+1:7*nd) = -25*d2r;
% ub1(4*nd+1:5*nd) = 1.04*omg_i*ones(nd,1);
% Gsol1P = Gsol1P(1:end-2);
lb1_new = [lb1(1:end-1); -8.4*d2r;  lb1(end)];
ub1_new = [ub1(1:end-1);  9.6*d2r; ub1(end)];
% lb1_new=lb1;ub1_new=ub1;
% lb1 = [lb1; -8.4*d2r; 1 - 1e-3];
% ub1 = [ub1;  8.4*d2r + 1e-3; 1];
% Gsol1P=Gsol1P(1:141);
if min(Gsol1P - lb1_new) <= 0
    lb1 = lb1 - 1e-3*ones(length(lb1),1);
end
if min(Gsol1P - ub1_new) >= 0
    ub1 = ub1 + 1e-3*ones(length(ub1),1);
end

% ub1(2*nd+1:3*nd) = ufmin + 1*0.5144;
% Gsol1P(nd+1:2*nd) = (Gsol1P(nd+1:2*nd)/1200)*800;
% Gsol1P(Gsol1P(4*nd+1:5*nd) > 1*omg_i) = 1*omg_i;
% Gsol1P(end) = 15;
D = Dmat(nd - 1);
[sol1P,~,flg1P] = fmincon(@costSSAR, Gsol1P, [],[],[],[], lb1_new, ub1_new, ...
                  @(x) NLC_thrPhs_P1(x, D, nd, BC1, heli, par), opts{1});
fprintf('Completed! Flag = %d\n', flg1P)
lb1(4*nd+1:5*nd) = 0.85*omg_i;
% ub1(2*nd+1:3*nd) = ufmin + 1*0.5144;
if -sol1P(2*nd) < 20/m2ft
    fprintf('Phase-1 Ended Under 20 ft.../n')
    s = input('');
    if isempty(s)
        return
    end
end

fprintf('Phase-2... ')
lb2 = [lb1(1:141); -8.4*d2r];
ub2 = [ub1(1:141);  9.6*d2r];
lb2(4*nd+1:5*nd) = 0.85*omg_i*ones(nd,1);

% BC2 = [sol1P([20,40]); ufmin; aTr(4); omg_i; aTr([3,2]); -20/m2ft];
BC2 = [sol1P(nd:nd:7*nd); -20/m2ft];
[sol2P,~,flg2P] = fmincon(@costLanding2, Gsol2P, [],[],[],[], lb2, ub2, ...
                  @(x) NLC_thrPhs_P2(x, D, nd, BC2, heli, par), opts{1});
fprintf('Completed! Flag = %d\n', flg2P)
fprintf('Phase-3... ')
lb3 = [lb1(1:120); -8.4*d2r*ones(nd,1); lb1(121:141)];
ub3 = [ub1(1:120);  9.6*d2r*ones(nd,1); ub1(121:141)];
lb3(4*nd+1:5*nd) = 0.60*omg_i*ones(nd,1); % Based on update on 300825

% Added Tail-Guard thing inside
% lb3(8*nd) = -7*d2r;
% ub3(8*nd) =  7*d2r;
BC3 = [sol2P(nd:nd:6*nd); sol2P(end); sol2P(7*nd); z_f];
[sol3P,~,flg3P] = fmincon(@costLanding, Gsol3P, [],[],[],[], lb3, ub3, ...
                  @(x) NLC_thrPhs_P3(x, D, nd, BC3, heli, par), opts{1});
fprintf('Completed! Flag = %d\n', flg3P)
fprintf('\nAll Flags: P1 = %1.0f P2 = %1.0f P3 = %1.0f\n', flg1P, flg2P, flg3P)

%%
nds   = chgslb(nd);
t3D   = 0.5*sol1P(end)*(nds + 1);
sol3D = [reshape(sol1P(1:6*nd), nd, 6), BC1(13)*ones(nd, 1), sol1P(6*nd+1:7*nd)];
sol3D = [sol3D; reshape(sol2P(1:6*nd), nd, 6), sol2P(end)*ones(nd, 1), sol2P(6*nd+1:7*nd)];
sol3D = [sol3D; reshape(sol3P(1:end-1), nd, 8)];

% Adding time-delay solution, if it was there
if ~all(par.TD == 0)
    t3D   = [tspan; t3D + tspan(end)];
    sol3D = [Xsol; sol3D];
end
t3D = [t3D; 0.5*sol2P(end-1)*(nds + 1) + t3D(end)];
t3D = [t3D; 0.5*sol3P(end)*(nds + 1) + t3D(end)];

% If the solution has to be saved for plotting/other purposes, set the
% variable "toSave" as 1. Otherwise, the code will stop run without saving.
% Change filename in variable "fname" if required.
toSave = 1;

if toSave == 1
    % Name of the .mat file
    fname = ['ThrPhs_MoI-',num2str(heli.rotor.MoI),'_v-',num2str(u_i),'_h-',num2str(z_i),'_TA(3.5,2).mat'];

    nTD = length(t3D) - 3*nd;
    save(fname, 'sol3D', 't3D', 'nd', 'nTD')
else
    return
end