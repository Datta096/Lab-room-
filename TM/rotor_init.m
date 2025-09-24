%% Conversions
ms2fpm  = 60/0.3048;   % m/s   -> ft/min
kts2ms  = 0.5144;      % kts   -> m/s
m2ft    = 1/0.3048;    % m     -> ft
rds2RPM = 30/pi;       % rad/s -> RPM
d2r     = pi/180;      % deg   -> rad

%% Earth Data
par.rho = 1.0879;      % Density
par.g   = 9.8;         % Acceleration due to gravity

%% Main Rotor Data
rotor.Nb     = 4;             % Number of rotor blades
rotor.R      = 5.8;           % Rotor disk radius
rotor.c      = 0.36;          % Rotor blade chord
cbar         = 0.36/rotor.R;  % Non-dimensional chord of the rotor
rotor.sigma  = 0.079;         % Rotor solidity
rotor.MB     = 45.1;          % Rotor Mass
rotor.e      = 1.425/rotor.R; % Rotor root cutout
rotor.tht_rt = 10.95*d2r;     % Rotor pitch at rbar = e
rotor.twist  = -0.1745;       % Rotor twist
rotor.Area   = pi*rotor.R^2;  % Rotor Disk Area
rotor.npsi   = 10;            % No. of azimuthal locations
rotor.alp_s  = -0.078525;     % Rotor shaft tilt
rotor.hngofs = 0.0841;        % Rotor blade hinge offset (e)
rotor.kb     = 5798;          % Rotor spring stiffness at the hinge
rotor.precon = 0.04363;       % Pre-cone Angle of the Main Rotor

% Aero Parameters
rotor.Cl0  = 0.0622; % Lift coefficient at zero angle of attack
rotor.Cla  = 5.8733; % Lift coefficient slope
rotor.Cd0  = 0.0058; % Drag coefficient at zero angle of attack
rotor.Cda  = 0.0229; % Drag coefficient linear term
rotor.Cda2 = 1.2732; % Drag coefficient square term

% Term for multiplying to coefficients
rotor.term = par.rho*rotor.Area*rotor.R^2;

% Legendre nodes for blade elements
nR             = 10; % No. of blade elements in a rotor blade
[xk, rotor.wk] = gauss(nR);
rotor.xk       = 0.5*((1 - rotor.e)*xk + (rotor.e + 1));

%% Other Control Surface Parameters
% H-Tail
ARH         = 5;
hvTRp.CD0H  = 0.079;
hvTRp.KH    = 1/(pi*0.8*ARH);
hvTRp.ClaH  = 5.73/(1 + 5.73*hvTRp.KH);
hvTRp.alpH0 = -0.02618;
hvTRp.SH    = 1.152;

% V-Tail
ARV         = 2.18;
hvTRp.CD0V  = 0.079;
hvTRp.KV    = 1/(pi*0.8*ARV);
hvTRp.ClaV  = 5.73/(1 + 5.73*hvTRp.KV);
hvTRp.alpV0 = 0.0698;
hvTRp.SV    = 0.96;

% Tail Rotor
hvTRp.sigTR   = 0.2405;
hvTRp.RTR     = 0.9;
hvTRp.TRGR    = 233.94/38.43;
hvTRp.del3    = 0.6494;
hvTRp.eTR     = 0.36/hvTRp.RTR;
hvTRp.thtRtTR = rotor.tht_rt + rotor.twist*(hvTRp.eTR - rotor.e);
hvTRp.twsTR   = rotor.twist;
hvTRp.gmaTR   = 1.79;
hvTRp.xk      = 0.5*((1 - hvTRp.eTR)*xk + (hvTRp.eTR + 1));

hvTRp.Cl0  = 0.0157;
hvTRp.Cla  = 4.8699;
hvTRp.Cd0  = 0.0070;
hvTRp.Cda  = 0.0130;
hvTRp.Cda2 = 0.3388;

%% Other Helicopter Data
% Dynamic Inflow Model Parameters
heli.M  = diag([128/(75*pi), 16/(45*pi), 16/(45*pi)]);
heli.t  = 15*pi/16;

% Inertia of the Helicopter
heli.Ix = 2255;
heli.Iy = 9756;
heli.Iz = 7870;

% CG-Shift
delCGx = 20*1e-2;
delCGy = -5*1e-2;

% Checking for variables... Ignore
if ~exist('delCGx','var'); delCGx = 0; end
if ~exist('delCGy','var'); delCGy = 0; end

% Locations
CG_SL = 4.194 + delCGx;
CG_WL = 2.299;
CG_BL = 0 + delCGy;

FS_SL = 4.162;
FS_WL = 2.299;
FS_BL = 0;

MR_SL = 4.162;
MR_WL = 3.902;
MR_BL = 0;

TR_SL = 10.891;
TR_WL = 2.918;
TR_BL = 0.385;

HF_SL = 9.378;
HF_WL = 2.7;
HF_BL = 0;

VF_SL = 11.424;
VF_WL = 3.3775;
VF_BL = 0;

heli.xFS = CG_SL - FS_SL;
heli.yFS = CG_BL - FS_BL;
heli.zFS = CG_WL - FS_WL;

heli.xMR = CG_SL - MR_SL;
heli.yMR = CG_BL - MR_BL;
heli.zMR = MR_WL - CG_WL;

heli.xHT = HF_SL - CG_SL;
heli.yHT = HF_BL - CG_BL;
heli.zHT = HF_WL - CG_WL;

heli.xVT = VF_SL - CG_SL;
heli.yVT = VF_BL - CG_BL;
heli.zVT = VF_WL - CG_WL;

heli.xTR = TR_SL - CG_SL;
heli.zTR = TR_WL - CG_WL;
heli.yTR = TR_BL - CG_BL;

% Fuselage Aerodynamics
load Aerofusdata.mat;
heli.CLfus = polyfit(CLfus(:,1), CLfus(:,2), 1);
heli.CDfus = polyfit(CLfus(:,1), CLfus(:,3), 2);
heli.CMfus = polyfit(CLfus(:,1), CLfus(:,4), 3);
heli.CMfus = heli.CMfus/4;
heli.CDfus = heli.CDfus*2;

%% Trajectory Optimization Related Parameters
nd = 20; % Number of CGL Nodes

% Boundary Conditions
t0    = 0.0;
w_i   = 0;
omg_i = 38.43;
x_i   = 0;
v_i   = 0;
u_i   = 0;
z_i   = -1200*0.3048;
z_f   = 0;
ufmin = 60*kts2ms;
 
%% Selection Parameters
Flag1 = 1; % Rotor Moment of Inertia
Flag2 = 1; % Thrust averaging: 0 - No, 1 - Yes
Flag3 = 0; % Non-linear model: 0 - No, 1 - Yes, 2 - sin()cos()
Flag4 = 0; % Time-delay      : 0 - No, 1 - Yes
Flag5 = 3; % Aircraft Mass

heli.f1 = Flag2;
heli.f2 = Flag3;
heli.f3 = 1; % Inflow Model: 1 -> New Fit, 0 -> Original Peters-He

% Rotor Moment of Inertia
if Flag1 == 1
    rotor.MoI  = 1550;
elseif Flag1 == 2
    rotor.MoI  = 1800;
elseif Flag1 == 3
    rotor.MoI  = 2500;
end

% Aircraft Mass
if Flag5 == 1
    heli.MF = 2500;
elseif Flag5 == 2
    heli.MF = 3150;
elseif Flag5 == 3
    % Custom
    heli.MF = 2400;
end

% Amount of time delay to account for pilot reaction time
if Flag4 == 0
    par.TD = 0;
else
    % Change value according to the requirement
    par.TD = 1;
end

% For collective dump. Set to the desired time (in sec). Otherwise, set to 0.
% Only use with main_3P_ind.m. Otherwise, it will lead to undesired effects.
CdumpTime = 1;
if CdumpTime ~= 0; par.TD = [par.TD, CdumpTime]; end

% Lock Number for the Main Rotor
rotor.gma  = par.rho*rotor.Cla*cbar*rotor.R^5/(rotor.MoI/rotor.Nb);

clear xk ARH ARV CG_SL CG_WL CG_BL FS_SL FS_WL MR_SL MR_WL TR_SL TR_WL...
    TR_BL HF_SL HF_WL VF_SL VF_WL CLfus Flag1 Flag2 Flag3 Flag4 Flag5 cbar delCGx...
    FS_BL MR_BL VF_BL HF_BL delCGy CdumpTime

%% Optimizer Bounds
ONS             = ones(nd, 1);
rotor.thtdotmax = deg2rad(8);
rotor.alpdotmax = deg2rad(15);
thtc_ofs        = 0*1;
thts_ofs        = 0*(-5);

% Lower bounds
xmin    =  -Inf*ONS;
zmin    = -Inf*ONS;
umin    =  0*ONS;
vmin    = -Inf*ONS;
wmin    = 0*ONS;
omgmin  = .80*omg_i*ONS;

thtmin  = -[25*ones(nd-1, 1); 25]*d2r;
thtDmin = -rotor.alpdotmax*ONS;
phiRmin = -[10*ones(nd-1, 1); 10]*d2r;
phiDmin = -rotor.alpdotmax*ONS;

lmb0min =    -0.2*ONS;
lmbsmin = -Inf*ONS;
lmbcmin = -Inf*ONS;
lmbTmin = -Inf*ONS;

tht0min = -15*d2r*ONS;
thtcmin = -(( 7 )*d2r)*ONS;
thtsmin = -((14 )*d2r)*ONS;
thtTmin = -(40*d2r)*ONS;

bet0min = -(20*d2r)*ONS;
betcmin = -(20*d2r)*ONS;
betsmin = -(20*d2r)*ONS;
betTmin = -(20*d2r)*ONS;

tfmin   = 0;

lb = [   xmin;    zmin;    umin;     wmin;
       omgmin; lmb0min;  thtmin;    tfmin];

% Upper bounds
xmax    = Inf*ONS;
zmax    = Inf*ONS;
umax    = Inf*ONS;
vmax    = Inf*ONS;
wmax    = Inf*ONS;
omgmax  = 1.05*omg_i*ONS;

thtmax  = [25*ones(nd-1, 1);25]*d2r;
thtDmax = rotor.alpdotmax*ONS;
phiRmax = [10*ones(nd-1, 1); 10]*d2r;
phiDmax = rotor.alpdotmax*ONS;

lmb0max = Inf*ONS;
lmbsmax = Inf*ONS;
lmbcmax = Inf*ONS;
lmbTmax = Inf*ONS;

tht0max = (15*d2r)*ONS;
thtcmax = ((9 - thtc_ofs)*d2r)*ONS;
thtsmax = ((7 - thts_ofs)*d2r)*ONS;
% thtcmax = (90*d2r)*ONS;
% thtsmax = (90*d2r)*ONS;
thtTmax = (40*d2r)*ONS;

bet0max = (20*d2r)*ONS;
betcmax = (20*d2r)*ONS;
betsmax = (20*d2r)*ONS;
betTmax = (20*d2r)*ONS;

tfmax   = 60;

ub = [   xmax;    zmax;    umax;    wmax;
       omgmax; lmb0max;  thtmax;   tfmax];

heli.rotor = rotor;
heli.hvTRp = hvTRp;

clear zmin wmin omgmin xmin umin vmin tht0min alfamin thtcmin thtsmin thtTmin...
      phiRmin bet0min betcmin betsmin betTmin lmb0min lmbsmin lmbcmin lmbTmin...
      zmax wmax omgmax xmax umax vmax tht0max alfamax thtcmax thtsmax thtTmax...
      phiRmax bet0max betcmax betsmax betTmax lmb0max lmbsmax lmbcmax lmbTmax...
      tfmin tfmax rotor hvTRp ONS phiDmax phiDmin thtDmax thtDmin thtmax thtmin...
      thtc_ofs thts_ofs