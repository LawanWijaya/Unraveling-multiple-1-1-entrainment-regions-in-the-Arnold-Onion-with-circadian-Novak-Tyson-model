% I am referring to the paper
% A Simple Model of Circadian Rhythms Based on Dimerization and
% Proteolysis of PER and TIM

% Phase Plane for squared Pulse
clear all;
close all;
clc
% Parameter values suitable for circadian rhythm of wild-type fruit flies
vm=1;
km=0.1;
vp=0.5;
kp1=10;
kp2=0.03;
kp3=0.1;
% Keq=200;
Pcrit=0.1;
Jp=0.05;
A=0.01;
php=0.5;
% per=60;

% setting values to the model runs
Keq=3.3;
% A=10;
per=63;

%time
tIni=0;
tend=per*5;
tspan=tIni:0.01:tend;

% ODE Solver
IC = [0.1 0.1];
options=odeset('AbsTol',1.e-8,'relTol',1.e-9,'InitialStep',1.e-3,'MaxStep',0.01);
[t1,y1]=ode15s(@NTdy,tspan,IC,options, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp, ...
    A,php,per);

IC = [y1(end,1) y1(end,2)];
[ti1,yy1]=ode15s(@NTdy,tspan,IC,options, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp, ...
    A,php,per);
% F = A*(1-heaviside(mod(ti1,per)-per*(1-php)));


M = yy1(:,1); %y(1)
Pt = yy1(:,2); %y(2)

% q = 2./(1+ sqrt(1+8.*Keq.*Pt));
% dMdt =@(P_t,M) vm./(1+(P_t.*(1-q)./(2.*Pcrit)).^2)-km.*M; %dM/dt
% dPtdt =@(P_t,M) vp.*M - (kp1.*P_t.*q + kp2.*P_t)./(Jp +P_t) - kp3.*P_t; %dPt/dt

%% For M nullcline when focing is 0
F=0; 
dMdt_0 =@(P_t,M) (1+F)*vm./(1+(P_t.*(1-(2./(1+ sqrt(1+8.*Keq.*P_t))))./(2.*Pcrit)).^2)-km.*M; %dM/dt

%% For M nullcline when focing is A
F=A; 
dMdt_A =@(P_t,M) (1+F)*vm./(1+(P_t.*(1-(2./(1+ sqrt(1+8.*Keq.*P_t))))./(2.*Pcrit)).^2)-km.*M; %dM/dt


dPtdt =@(P_t,M) vp.*M - (kp1.*P_t.*(2./(1+ sqrt(1+8.*Keq.*P_t))) + kp2.*P_t)./(Jp +P_t) - kp3.*P_t; %dPt/dt


%% Give red or black to limit cycle based on F =A or F=0
% F_values = A*(1-heaviside(mod(ti1,per)-per*(1-php))); % php value shows "off" percentage
F_values = A*(heaviside(mod(ti1,per)-per*(1-php))); % php value shows "on" percentage
% F_values = A*(heaviside(per*php-mod(ti1,per))); % similar outcome to the above line

% Initialize arrays for different forcing states
M_forcing_on = M; M_forcing_on(F_values==0) = nan;  % Red segment
M_forcing_off = M; M_forcing_off(F_values>0) = nan; % Black segment
Pt_forcing_on = Pt; Pt_forcing_on(F_values==0) = nan;
Pt_forcing_off = Pt; Pt_forcing_off(F_values>0) = nan;

% input and the response
figure(1)
subplot(3,1,[1 2])
% plot(ti1,M,LineWidth=2)
% Plot M in black when F_values is 0
plot(ti1, M_forcing_off, 'k', 'LineWidth', 2.5)
hold on
% Plot M in red when F_values is non-zero
plot(ti1, M_forcing_on, 'r', 'LineWidth', 2.5)
hold on
% plot(ti1,Pt,Color=[1 0 0],LineWidth=2)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',21)
% xlabel('Time','FontSize',18,'Interpreter','latex')
ylabel('$M$','FontSize',21,'Interpreter','latex')
% legend('mRNA','Protein','interpreter','latex')
% xlim([0 80]);
xlim([0 tend]);


subplot(3,1,3)
plot(ti1,F_values,'-r','linewidth',2.5)
hold on
% Overlay the zero parts with black
zeroSegments = (F_values == 0);


% Plot the zero segments in black
% Use the find function to identify start and end of zero segments
startIdx = find(diff(zeroSegments) == 1) + 1; % Start of zero segments
endIdx = find(diff(zeroSegments) == -1); % End of zero segments

% For each zero segment, plot a short black line
for k = 1:length(startIdx)
    % Ensure indices are within bounds
    if k <= length(endIdx)
        plot(ti1(startIdx(k):endIdx(k)), F_values(startIdx(k):endIdx(k)), 'k', 'LineWidth', 2);
    end
end

% Adjust the plot settings
% hold off;
xlabel('Time', 'Interpreter', 'latex');
ylabel('F', 'Interpreter', 'latex');
xlim([0 tend]);
% ylim([min(F_values) - 0.0001, max(F_values) + 0.006]);
ylim([-0.01 0.41])
set(gca, 'FontSize', 21);

% Plot F_values in black when it's zero
% plot(ti1(F_values == 0), F_values(F_values == 0), 'k', 'LineWidth', 2)
% hold on
% % Plot F_values in red when it's non-zero
% plot(ti1(F_values > 0), F_values(F_values > 0), 'r', 'LineWidth', 2)



set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'FontSize',21)
yticks([0, max(F_values)]);
xlabel('Time','Interpreter','latex')
ylabel('$F$','Interpreter','latex')
xlim([0 tend]);
ylim([min(F_values)-0.0001 max(F_values)+0.006])
box
hold on
box

% phase plane
figure(2)
hold on;

% fimplicit(dMdt,[0.001 100,0.1 10],'LineWidth',2,'Color','r');
fimplicit(dMdt_A,[0.001 100,1 15],'LineWidth',2.5,'Color','r','LineStyle','--');
hold on;
fimplicit(dMdt_0,[0.001 100,1 15],'LineWidth',2.5,'Color','k','LineStyle','--');
fimplicit(dPtdt,[0.001 100,1 15],'LineWidth',2.5,'Color','g');
hold on;
% plot(Pt,M,'LineWidth',2,'Color','r')
% Plot Forcing On in Red
plot(Pt_forcing_on, M_forcing_on, 'r', 'LineWidth', 2.5);
% Plot Forcing Off in Black
plot(Pt_forcing_off, M_forcing_off, 'k', 'LineWidth', 2.5);
% hold off;

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',21)
% xlabel('Protein','FontSize',18,'Interpreter','latex')
xlabel('$P_t$','FontSize',21,'Interpreter','latex')
% ylabel('mRNA','FontSize',18,'Interpreter','latex')
ylabel('$M$','FontSize',21,'Interpreter','latex')
% % % title(['A = ', num2str(A), ', T = ', num2str(per)], 'FontSize', 18, 'Interpreter', 'latex');
% legend('M nullcline (F = A)','M nullcline (F = 0)','$P_t$ nullcline','Limit cycle','interpreter','latex')
xlim([0.001 100]);
box


function dy = NTdy(ti1,y, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp, ...
    A,php,per)
F = A*(1-heaviside(mod(ti1,per)-per*(1-php))); % square pulse
q = 2/(1+ sqrt(1+8*Keq*y(2)));
% M = y(1)
% Pt = y(2)
dy(1,1) = (1+F)*vm/(1+(y(2)*(1-q)/(2*Pcrit))^2)-km*y(1); %dM/dt
dy(2,1) = vp*y(1) - (kp1*y(2)*q + kp2*y(2))/(Jp +y(2)) - kp3*y(2); %dPt/dt

end