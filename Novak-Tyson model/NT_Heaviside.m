% Manual Arnold tongue finding
% Using mRNA for spikes

% I am referring to the Thsise
% A Simple Model of Circadian Rhythms Based on Dimerization and
% Proteolysis of PER and TIM

clear all;
close all;
clc

%% Parameter values suitable for circadian rhythm of wild-type fruit flies
vm=1; % this is vm in the autonomous one % maximum rate of synthesis of mRNA
km=0.1; % First-order rate constant for mRNA degradation
vp=0.5; % Rate constant for translation of mRNA
kp1=10; % Vmax for monomer phosphorylation
kp2=0.03; % Vmax for dimer phosphorylation
kp3=0.1; % First-order rate constant for proteolysis
% Keq=200; % Equilibrium constant for dimerization

Pcrit=0.1; % Dimer concen at the half-maximum transcription rate
Jp=0.05; % Michaelis constant for protein kinase (DBT)


%% For time dependent squared pulse forcing and changing variables
A=2;
% vm=1; % 0.2208 and 2.858 are bifurcation points
Keq=2; % equilibrium constant for dimerization
% per=17; % no impact when A =0
% php=0.5; % no impact when A=0

%time

tIni=0;
php=0.5;
per=24; 
tstart=per*0;
% tend=250.8;
tend=per*5;
% tend=50;
tspan=tIni:0.01:tend;

%% ODE Solver
IC = [0.1 0.1];
options=odeset('AbsTol',1.e-5,'relTol',1.e-5,'InitialStep',1.e-3,'MaxStep',0.01);
[t1,y1]=ode15s(@NTdy,tspan,IC,options, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per);
IC = [y1(end,1) y1(end,2)];
[ti1,yy1]=ode15s(@NTdy,tspan,IC,options, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per);


% IC = [y1(end,1) y1(end,2)];
% [ti1,yy1]=ode15s(@NTdy,tspan,IC,options, ...
%     vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per);


F = A*(heaviside(php*per-mod(ti1,per))); % square pulse
% F=A*(1+sign(sin((2*pi*ti1)/per)))/2;

M = yy1(:,1); %y(1)
Pt = yy1(:,2); %y(2)
%% counting
% Find peaks in mRNA levels
[pks, locs] = findpeaks(M,'MinPeakProminence',(max(M)-min(M))/2);

% Filter peaks and pulses within the time interval from 300 to 600
valid_locs = locs(locs >= tstart*100);% & locs <= 600);
valid_pks = pks(locs >= tstart*100);% & locs <= 600);

% Count the number of spikes within the interval
num_spikes = numel(valid_locs);

% Find the number of pulses (periodic peaks) within the interval
% mean_period = mean(diff(valid_locs));
% mean_period=(tend-tstart)/(per*(1+php));
num_pulses = floor((tend - tstart) / per);
ratio=num_spikes/num_pulses;
disp(['Number of spikes within the interval: ', num2str(num_spikes)]);
disp(['Number of pulses within the interval: ', num2str(num_pulses)]);
disp(['Ratio: ', num2str(ratio)]);

% disp(['Ratio'])
%%
figure(1)
subplot(3,1,[1 2])
plot(ti1-tstart,M,LineWidth=2)
% hold on
% plot(ti1,Pt,Color=[1 0 0])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'FontSize',24)
% xlabel('Time','Interpreter','latex')
% ylabel('[mRNA]','Interpreter','latex')
ylabel('$M$','Interpreter','latex')
% legend('mRNA','Protein','interpreter','latex')
xlim([0 tend-tstart]);
% ylim([min(M)-1 max(M)+2])
ylim([0 max(M)+2]);
set(gca, 'XTickLabel', []);

subplot(3,1,3)
plot(ti1-tstart,F,'-r','linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'FontSize',24)
yticks([0, max(F)]);
xlabel('Time','Interpreter','latex')
ylabel('$F$','Interpreter','latex')
xlim([0 tend-tstart]);
ylim([min(F)-0.001 max(F)+A])

% legend('mRNA','Protein','interpreter','latex')
% xlim([0 tend-tstart]);
% ylim([vm-0.0001 vm+A*0.1+0.002])

% % % figure(2)
% % % plot(ti1-tstart,F,'-r','linewidth',2)
% % % set(groot,'defaultAxesTickLabelInterpreter','latex');
% % % set(groot,'defaulttextinterpreter','latex');
% % % set(groot,'defaultLegendInterpreter','latex');
% % % set(gca,'FontSize',24)
% % % yticks([0, max(F)]);
% % % xlabel('Time','Interpreter','latex')
% % % ylabel('$F$','Interpreter','latex')
% % % xlim([0 tend-tstart]);
% % % ylim([min(F)-0.5 max(F)+1])

function dy = NTdy(ti1,y,vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per)

F = A*(1-heaviside(mod(ti1,per)-per*(1-php))); % square pulse
% F=A*(1+sign(sin((2*pi*ti1)/per)))/2;

% M = y(1)
% Pt = y(2)
q = 2/(1+ sqrt(1+8*Keq*y(2)));
dy(1,1) = (1+F)*vm/(1+(y(2)*(1-q)/(2*Pcrit))^2)-km*y(1); %dM/dt
dy(2,1) = vp*y(1) - (kp1*y(2)*q + kp2*y(2))/(Jp +y(2)) - kp3*y(2); %dPt/dt

end