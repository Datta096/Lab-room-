clc; clear; close all;

rotor_init;
clr = {'r'; 'b'; 'k'; 'm'; 'g'};
fil = repmat({'-','--',':','-.'}, length(clr), 1);
cfl = [repmat(clr, size(fil,2), 1), fil(:)];
clear clr fil

% Loading the data
[lgd, folderPath] = uigetfile('*.mat','Select mat files','Multiselect','on');
%disp(lgd)
if ischar(lgd) 
    lgd = cellstr(lgd);
end

for k = 1:length(lgd)
    stemp = lgd(k);
    S{k} = load(stemp{1});
    
    % Post-processing inputs' solutions
    S{k}.sol6DE(:,15) = (S{k}.sol6DE(:,15)/d2r + 8.4)*100/18; % Collective  (theta0)
    S{k}.sol6DE(:,16) = (S{k}.sol6DE(:,16)/d2r + 7.0)*100/16; % Lat. Cyclic (thetac)
    S{k}.sol6DE(:,17) = (7.0 - S{k}.sol6DE(:,17)/d2r)*100/21; % Lng. Cyclic (thetas)
end

nS = length(S);


% If the solution had time-delay, then the solution saved in the file would
% have the time vector starting from par.TD. If this is not required, then
% uncomment the following line.
% S.t6DE = S.t6DE - S.t6DE(1);

% Figure width and height in pixels
Pix_SS = get(0,'screensize');
width  = 1200; height = 900;
fsize  = [(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height];



mult = [m2ft, -m2ft, 1/kts2ms, ms2fpm, 100/38.43, 1/d2r, 1, 1, 1];
tlt  = {'x (in ft)', 'z (in ft)', 'u (in kts)', 'w (in fpm)', '\Omega (in %)',...
        '\theta_P (in deg)', '\theta_0 (in %)', '\theta_c (in %)', '\theta_s (in %)'};
ind  = [1:3, 5:7, 15:17];

%%
close all;

f1 = figure(1);
f1.Position = fsize;
for j = 1:nS
    Si = S{j};
    cl = cfl{mod(j,length(cfl)), 1};
    fl = cfl{mod(j,length(cfl)), 2};
    for k = 1:9
        subplot(3,3,k)
        plot(Si.t6DE, Si.sol6DE(:,ind(k))*mult(k), 'Color', cl, 'LineStyle', fl, 'LineWidth', 1,'marker', '*')
        grid on
        set(gca,'FontName','Lucida Sans','FontSize',12)
        xlabel('Time (in sec)', 'FontSize',14)
        title(tlt{k}, 'FontSize',14)
        if j==1;hold on;end 
        
    end
    
end
subplot(3,3,1)
legend(lgd, 'FontSize',14, 'Location','best','interpreter','none')