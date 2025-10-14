clc; clear; close all;
rotor_init;
par.p=0;
opts   = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12, ...
                      'Algorithm','trust-region-dogleg','Display','none', ...
                      'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);

% adding xls read part
    [xlfile, folderPath] = uigetfile('*.xlsx','Select xlsx files', 'D:\LUH\PT2\improved inertia MR\H-V related\Takeoff corridor\' ,'Multiselect','off'); % kindly change the default folder if any error occurs
    pdata = [ xlsread([folderPath xlfile],'P','AH3:AH783') xlsread([folderPath xlfile],'P','t3:t783') xlsread([folderPath xlfile],'P','Aj3:Aj783') xlsread([folderPath xlfile],'P','Ae3:Ae783') xlsread([folderPath xlfile],'P','B3:B783')];
    pdata = [pdata xlsread([folderPath xlfile],'P','D3:D783') xlsread([folderPath xlfile],'P','C3:C783') xlsread([folderPath xlfile],'P','V3:V783')]; % Adding cyclic also
    Ddata = [xlsread([folderPath xlfile],'D','J3:J302') xlsread([folderPath xlfile],'D','M3:P302')];
    refht =  xlsread([folderPath xlfile],'D','t1')*3.28084;
% end of xls read part
% load Pdata1.mat % Time Pitch Torque Omega Collective %thtc thts% roll (New)
% load Ddata1.mat % X t U W Z
data1=pdata;
 %data1(:,3)=Pdata(:,3)*Tq_max/100;
for i=1:8 %5
    data2(:,i)=interp1(data1(:,1),data1(:,i),Ddata(:,2));
end
xdata=data2;
data2(:,5)=((data2(:,5)*18/100)-8.4)*d2r;
data2(:,2)=data2(:,2)*d2r;
data2(:,4)=data2(:,4)*omg_i/100;

data2(:,6)=((data2(:,6)*16/100)-7)*d2r;
data2(:,7)=(-(data2(:,7)*21/100)+7)*d2r;

data3=data2;
data3=[data3(:,1:end-1),Ddata(:,1),Ddata(:,3:end),data3(:,end)*d2r];

% Changing w into sink rate
% data3(:,7)=-data3(:,7);
data3(:,10) = -data3(:,10);

% data3(:,3) = zeros(300,1); % as a way to reduce the torque to zero

% cutting Data3 variable to start at enG OFF and end at touchdown
ENGoff = xlsread([folderPath xlfile],'Results','F7'); % getting eng off time from xls file
Touchdown = xlsread([folderPath xlfile],'Results','G7'); % getting touchdown time from xls file
EV = round(ENGoff,1)*10; % rounding to 1 decimal place so as to match with data3 first column
TD = round(Touchdown,1)*10;
data3 = data3(EV:TD,:); % reducing data3 to only ENG off to touchdown time
data3(:,1) = data3(:,1) - data3(1,1); % adjusting time column to start from zero.
%data3 Time Pitch Torque Omega Collective %thtc thts% (New) X U W Z roll
data4=[0.5*data3(end,1)*(chgslb(20) + 1)];

for i = 2:12
    data4 = [data4, interp1(data3(:,1), data3(:,i), data4(:,1))];
end

% for i=1:40
%     if rem(i,2)~=0
%         data4=[data4;data3(i,:)];
%     end
% end
% Initial x-location (in ft)
x_i = 0;

% Initial z-location (altitude, in ft)
z_i = data4(1,11)-refht; % changed ground altitude from 2935 to reference height given in xlsfile

% Initial forward velocity (in kts)
u_i = data4(1,9);

% Initial sink rate (in fpm)
w_i = data4(1,10);

% Initial rotor angular speed (in %)
o_i = data4(1,4)*100/38.43;
sol3D=[data4(:,8),-(data4(:,11)-refht)/m2ft,data4(:,9)*kts2ms,data4(:,10)/ms2fpm,data4(:,4)];

function out = tmpfunc( sts, ctrl,lmb, heli, par)
    dx = LDM_state_eqns(sts, [ctrl; lmb], heli, par, 2);
    out = dx(3);
end

for i=1:20
    sts = sol3D(i,:);
    ctrl=[data4(i,5);data4(i,2)];
    [lam(i),~,flg] = fsolve(@(lmb) tmpfunc(sts, ctrl,lmb, heli, par), 0.01, opts);
    if ~any(flg == 1:4); fprintf('lambda_i not solved i = %d\n', i); end
end
S.t3D=data4(:,1);
S.sol3D=[data4(:,8),-(data4(:,11)-refht)/m2ft,data4(:,9)*kts2ms,data4(:,10)/ms2fpm,data4(:,4),lam',data4(:,5),data4(:,2),data4(:,end)];
save fltdata.mat data4
S.nd=20;S.nTD=0;
% Initializing Parameters

gusfsl = [0.1; 0.1; 0.001];
sts = [x_i/m2ft; -z_i/m2ft; u_i*kts2ms; w_i/ms2fpm; o_i*38.43/100];

[out,~,flg] = fsolve(@(ctrl) LDM_state_eqns(sts, ctrl, heli, par, 2), gusfsl,opts );

if any(flg == 1:4)
    fprintf('Trim Solution Obtained...\n')
end
tht0 = out(1);
thtP = out(2);
l_i  = out(3);
sts1=[sts;l_i];
ctrl1=[tht0;thtP];
x1=LDM_state_eqns(sts1,ctrl1,heli,par,1);

% Optimizer Parameters
opts = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12, ...
                    'Algorithm','trust-region-dogleg','Display','none', ...
                    'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);


% Finding Autorotation Trim of 6-DoF system to use as a guess
load trim_AR_guess_20.mat
gusart = out;

trARF = @(x, sts) state_eqns_6D(sts(1:4), [x; sts(5)], heli, par, 4);

[aTr,~,flg] = fsolve(@(x) trARF(x, sts), gusart, opts);

if any(flg == [1,2,3,4])
    fprintf('Autorotation Trim Solution Obtained\n')
end

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
    attd(:,4) = sol(:,end);   % phiR
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
    fname = 'sol6DE_fd.mat';

    save(fname, 'sol6DE', 't6DE', 'nd', 'nP')
else
    return
end
