clc; clear; close all;

rotor_init;
clr = {'r'; 'b'; 'k'; 'm'; 'g'};
fil = repmat({'-','--',':','-.'}, length(clr), 1);
cfl = [repmat(clr, size(fil,2), 1), fil(:)];
clear clr fil

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% README %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This code is for plotting and comparing various individual run solutions.
Each solution is taken from a separate file.

All the loaded solution should have the same variable names, viz., sol3D and t3D.
1. sol3D should be a matrix with 8 columns [x, z, u, w, omg, lambda0, theta0, thetaP]
2. t3D should be a vector, with the same length as sol3D.

Please pre-process the data to ensure this is the case. Otherwise, errors will arise.
All the code written in other main files will save it in this format. So, no problem.

Otherwise, copy-paste the code for plotting elsewhere and use it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Load all the files here. Save it as S{number}, number being a value based
% on the order of loading. First file is S{1} and so on.
%S{1} = load('ThrPhs_v-0_h-1200.mat');
% Fill legends here in the order of data
%lgd = {'ThrPhs_v-0_h-1200_AUW-2800_MOI-1550.mat','ThrPhs_v-0_h-1200_AUW-2800_MOI-1800.mat','ThrPhs_v-0_h-1200_AUW-2800_MOI-2000.mat','ThrPhs_v-0_h-1200_AUW-2800_MOI-2500.mat'};
%lgd = {'ThrPhs_v-0_h-1200_AUW-2800_MOI-1800.mat','ThrPhs_v-0_h-1200_AUW-2800_MOI-1800_Radius-5.7_chord-0.34.mat','ThrPhs_v-0_h-1200_AUW-2800_MOI-1800_Radius-5.6_chord-0.34.mat'}
%lgd = {'ThrPhs_v-0_h-2000_AUW-2800_MOI-2500.mat'}

[lgd, folderPath] = uigetfile('*.mat','Select mat files','Multiselect','on');
%disp(lgd)
if ischar(lgd) 
    lgd = cellstr(lgd);
end

for i = 1:length(lgd)
    stemp = lgd(i);
    S{i} = load(stemp{1});
end

nS = length(S);

% Some post-processing
for i = 1:nS
    S{i}.sol3D(:,7) = (S{i}.sol3D(:,7)/d2r + 8.4)*100/18;
end

% Figure width and height in pixels
Pix_SS = get(0,'screensize');
width  = 1200; height = 900;
fsize  = [(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height];

mult = [m2ft, -m2ft, 1/kts2ms, ms2fpm, 100/38.43, 1, 1, 1/d2r];
tlt  = {'x (in ft)', 'z (in ft)', 'u (in kts)', 'w (in fpm)', '\Omega (in %)',...
        '\lambda_i', '\theta_0 (in %)', '\theta_P (in deg)'};

%%
close all;
f1 = figure(1);
f1.Position = fsize;

for i = 1:nS
    Si = S{i};
    cl = cfl{mod(i,length(cfl)), 1};
    fl = cfl{mod(i,length(cfl)), 2};
    for j = 1:8
        subplot(3,3,j)
        plot(Si.t3D, Si.sol3D(:,j)*mult(j), 'Color', cl, 'LineStyle', fl, 'LineWidth', 1)
        if i == 1; hold on; end
    end
end

for i = 1:8
    subplot(3,3,i)
    grid on
    set(gca,'FontName','Lucida Sans','FontSize',12)
    xlabel('Time (in sec)','FontSize',14)
    title(tlt{i},'FontSize', 14)
end

subplot(3,3,9)
for i = 1:nS
    cl = cfl{mod(i,length(cfl)), 1};
    fl = cfl{mod(i,length(cfl)), 2};
    plot(nan, nan, 'Color', cl, 'LineStyle', fl, 'LineWidth', 1)
    if i == 1; hold on; end
   
end
legend(lgd, 'FontSize',14, 'Location','best','interpreter','none')
set(gca,'XColor', 'none','YColor','none')

%%
%{
% The following commented code uses tiledlayout. If required, uncomment by
% removing the curly bracket { from the line above.

close all;
f1 = figure(1);
f1.Position = fsize;
t = tiledlayout(3,3, 'TileSpacing','tight', 'Padding','tight');

for i = 1:nS
    Si = S{i};
    cl = cfl{mod(i,length(cfl)), 1};
    fl = cfl{mod(i,length(cfl)), 2};
    for j = 1:8
        nexttile(j)
        plot(Si.t3D, Si.sol3D(:,j)*mult(j), 'Color', cl, 'LineStyle', fl, 'LineWidth', 1)
        if i == 1; hold on; end
    end
end

for i = 1:8
    nexttile(i)
    grid on
    set(gca,'FontName','Lucida Sans','FontSize',12)
    title(tlt{i},'FontSize', 14)
end

xlabel(t,'Time (in sec)', 'FontSize',16)
lgd = legend(lgd{1}, 'FontSize', 14);
lgd.Layout.Tile = 9;

%}