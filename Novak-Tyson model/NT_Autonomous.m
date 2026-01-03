% I am referring to the paper
% A Simple Model of Circadian Rhythms Based on Dimerization and
% Proteolysis of PER and TIM

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
Keq=200;
Pcrit=0.1;
Jp=0.05;

%time
tIni=0;
tend=100;
tspan=tIni:0.01:tend;

% ODE Solver
IC = [0.1 0.1];
options=odeset('AbsTol',1.e-8,'relTol',1.e-9,'InitialStep',1.e-3,'MaxStep',0.01);
[t1,y1]=ode15s(@NTdy,tspan,IC,options, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp);

IC = [y1(end,1) y1(end,2)];
[ti1,yy1]=ode15s(@NTdy,tspan,IC,options, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp);

M = yy1(:,1); %y(1)
Pt = yy1(:,2); %y(2)

% q = 2./(1+ sqrt(1+8.*Keq.*Pt));
% dMdt =@(P_t,M) vm./(1+(P_t.*(1-q)./(2.*Pcrit)).^2)-km.*M; %dM/dt
% dPtdt =@(P_t,M) vp.*M - (kp1.*P_t.*q + kp2.*P_t)./(Jp +P_t) - kp3.*P_t; %dPt/dt
dMdt =@(P_t,M) vm./(1+(P_t.*(1-(2./(1+ sqrt(1+8.*Keq.*P_t))))./(2.*Pcrit)).^2)-km.*M; %dM/dt
dPtdt =@(P_t,M) vp.*M - (kp1.*P_t.*(2./(1+ sqrt(1+8.*Keq.*P_t))) + kp2.*P_t)./(Jp +P_t) - kp3.*P_t; %dPt/dt


figure(1)
plot(ti1,M)
hold on
plot(ti1,Pt,Color=[1 0 0])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18)
xlabel('Time (hr)','FontSize',18,'Interpreter','latex')
ylabel('mRNA \& Protein Level','FontSize',18,'Interpreter','latex')
legend('mRNA','Protein','interpreter','latex')
xlim([0 80]);


figure(2)
fimplicit(dMdt,[0.001 100,0.1 10],'LineWidth',2,'Color','r');
hold on;
fimplicit(dPtdt,[0.001 100,0.1 10],'LineWidth',2,'Color','b');
hold on;
plot(Pt,M,'LineWidth',2,'Color','g')
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18)
xlabel('Protein','FontSize',18,'Interpreter','latex')
ylabel('mRNA','FontSize',18,'Interpreter','latex')
legend('M nullcline','$P_t$ nullcline','Limit cycle','interpreter','latex')

figure(3)
fimplicit(dMdt,[0.001 100,0.1 10],'LineWidth',2,'Color','r');
hold on;
fimplicit(dPtdt,[0.001 100,0.1 10],'LineWidth',2,'Color','b');
hold on;
plot(Pt,M,'LineWidth',2,'Color','g')
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18)
xlabel('Protein','FontSize',18,'Interpreter','latex')
ylabel('mRNA','FontSize',18,'Interpreter','latex')
legend('M nullcline','$P_t$ nullcline','Limit cycle','interpreter','latex')

% Now let's include the direction field

x_values = logspace(-2, 2, 20);
y_values = logspace(-1, 1, 20);
[x, y] = meshgrid(x_values, y_values); % Define grid for plotting vectors
u = zeros(size(x));
v = zeros(size(y));
for i = 1:numel(x)
    dydt = NTdy(0, [x(i); y(i)], vm, km, vp, kp1, kp2, kp3, Keq, Pcrit, Jp);
    u(i) = dydt(1);
    v(i) = dydt(2);
end
quiver(x, y, u, v, 'Color', 'k', 'LineWidth', 1.5);

function dy = NTdy(t,y, ...
    vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp)

q = 2/(1+ sqrt(1+8*Keq*y(2)));
% M = y(1)
% Pt = y(2)
dy(1,1) = vm/(1+(y(2)*(1-q)/(2*Pcrit))^2)-km*y(1); %dM/dt
dy(2,1) = vp*y(1) - (kp1*y(2)*q + kp2*y(2))/(Jp +y(2)) - kp3*y(2); %dPt/dt

end