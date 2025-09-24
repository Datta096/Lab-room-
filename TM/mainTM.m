clc; clear; close all;
rotor_init;

opts   = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12, ...
                      'Algorithm','trust-region-dogleg','Display','none', ...
                      'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);

load Pdata.mat
load Ddata.mat
data1=Pdata;
% data1(:,3)=Pdata(:,3)*Tq_max/100;
for i=1:5
    data2(:,i)=interp1(data1(:,1),data1(:,i),Ddata(:,1));
end
xdata=data2;
data2(:,5)=((data2(:,5)*18/100)-8.4)*d2r;
data2(:,2)=data2(:,2)*d2r;
data2(:,4)=data2(:,4)*omg_i/100;
data3=data2;
data3=[data3,Ddata(:,2:end)];
data3(:,7)=-data3(:,7);

ind = find((data3(:,1) >= 1.04) & data3(:,1) <= 5.65);
data3 = data3(ind,:);
xdata = xdata(ind,:);
data3(:,1) = data3(:,1) - data3(1,1);
xdata(:,1) = xdata(:,1) - xdata(1,1);

% Initial x-location (in ft)
x_i = 0;

% Initial z-location (altitude, in ft)
z_i = data3(1,8);

% Initial forward velocity (in kts)
u_i = data3(1,6);

% Initial sink rate (in fpm)
w_i = data3(1,7);

% Initial rotor angular speed (in %)
o_i = data3(1,4)*100/38.43;

sts = [x_i/m2ft; -z_i/m2ft; u_i*kts2ms; w_i/ms2fpm; o_i*38.43/100];
% Initializing Parameters

gusfsl = [0.1; 0.1; 0.001];

[out,~,flg] = fsolve(@(ctrl) LDM_state_eqns(sts, ctrl, heli, par, 2,0), gusfsl,opts );

if any(flg == 1:4)
    fprintf('Trim Solution Obtained...\n')
end
tht0 = out(1);
thtP = out(2);
l_i  = out(3);
sts1=[sts;l_i];
ctrl1=[tht0;thtP];
x1=LDM_state_eqns(sts1,ctrl1,heli,par,1,0);
Tq_max=-x1(5)*heli.rotor.MoI;
% Tq_max=0;
data3(:,3)=data3(:,3)*Tq_max/100;

P1={data3(:,1),data3(:,5)};
P2={data3(:,1),data3(:,2)};
P3={data3(:,1),data3(:,3)};
P4={data3(:,1),data3(:,4)};
tsp=data3(:,1);
X0=sts1;
odeopt = odeset('RelTol',1e-10, 'AbsTol',1e-10,'Events',@(t,x)myE(t,x,P2));
function [v, t, d] = myE(t, x, P2)%, lnsP, tmsP)
    thtP = interp1(P2{1}, P2{2}, t);
    v = x(2) + (74.956*thtP*180/pi + 1502.6)*1e-3;
    t = 1;
    d = 1;
end
[tsp1, Xsol,te,~] = ode45(@(t,x) tmFunc(t,x, P1, P2,P3, heli, par), tsp, X0([1:4,6]), odeopt);

% ind = find(data3(:,1) <= tsp1(end));
% data3 = data3(ind,:);
% xdata = xdata(ind,:);
% Xsol=[Xsol(:,1:4),zeros(length(tsp1),1),Xsol(:,5)];
figure(1)
tiledlayout(3,2)
nexttile
plot(tsp1, -Xsol(:,2)*m2ft, 'r*-')
hold on
plot(data3(:,1),data3(:,8),'b*-')
grid on
set(gca,'FontName','Lucida Sans','FontSize',12)
xlabel('Time (in sec)', 'FontSize',14)
title('z (in ft)', 'FontSize',14)
nexttile
plot(tsp1, Xsol(:,3)/kts2ms, 'r*-')
hold on
plot(data3(:,1),data3(:,6),'b*-')
grid on
set(gca,'FontName','Lucida Sans','FontSize',12)
xlabel('Time (in sec)', 'FontSize',14)
title('u (in kts)', 'FontSize',14)
nexttile
plot(tsp1, Xsol(:,4)*ms2fpm, 'r*-')
hold on
plot(data3(:,1),data3(:,7),'b*-')
grid on
set(gca,'FontName','Lucida Sans','FontSize',12)
xlabel('Time (in sec)', 'FontSize',14)
title('w (in fpm)', 'FontSize',14)
nexttile
plot(tsp1, Xsol(:,5)*100/omg_i, 'r*-')
hold on
plot(data3(:,1),xdata(:,4),'b*-')
grid on
set(gca,'FontName','Lucida Sans','FontSize',12)
xlabel('Time (in sec)', 'FontSize',14)
title('omega (in %)', 'FontSize',14)
nexttile
plot(data3(:,1),(data3(:,5)*180/pi+8.4)*100/18,'b*-')
grid on
set(gca,'FontName','Lucida Sans','FontSize',12)
xlabel('Time (in sec)', 'FontSize',14)
title('thet0 (in %)', 'FontSize',14)
nexttile
plot(data3(:,1),data3(:,2)*180/pi,'b*-')
grid on
set(gca,'FontName','Lucida Sans','FontSize',12)
xlabel('Time (in sec)', 'FontSize',14)
title('thetP (in deg)', 'FontSize',14)
