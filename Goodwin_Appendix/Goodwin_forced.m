% Goodwin Oscillator - With square forcing
clearvars;
close all;
clc

per=3.6; %0 period of the forcin7
% per=4.5; %0 period of the forcin7
Total_time=per*30;
T_shift=per*10;
time_step=0.01;
time=0:time_step:Total_time;

% Initilaization
x=zeros(1,length(time));
y=zeros(1,length(time));
z=zeros(1,length(time));
f=zeros(1,length(time));

% Iinitial values of state variables
x(1)=0.1;
y(1)=0.1;
z(1)=0.1;


alpha=100;


php=0.5; % photoperiod (or duty cycle)

n=7.5;

Amp=0.1; % strength of the forcing


tm=mod(time,per);
tmtest=tm-per*(1-php);
Forcing=Amp*heaviside(tmtest); % square pulses
for j=1:length(time)-1


    x1=alpha*(1+Forcing(j))*(1/(1+z(j)^n))-x(j);
    y1=x(j)-y(j);
    z1=y(j)-z(j);

    incx=x(j)+x1*time_step;
    incy=y(j)+y1*time_step;
    incz=z(j)+z1*time_step;


    x2=alpha*(1+Forcing(j+1))*(1/(1+incz^n))-incx;
    y2=incx-incy;
    z2=incy-incz;

    x(j+1)=x(j)+(x1+x2)*time_step/2;
    y(j+1)=y(j)+(y1+y2)*time_step/2;
    z(j+1)=z(j)+(z1+z2)*time_step/2;
end

x(1)=x(end);
y(1)=y(end);
z(1)=z(end);
for j=1:length(time)-1


    x1=alpha*(1+Forcing(j))*(1/(1+z(j)^n))-x(j);
    y1=x(j)-y(j);
    z1=y(j)-z(j);

    incx=x(j)+x1*time_step;
    incy=y(j)+y1*time_step;
    incz=z(j)+z1*time_step;


    x2=alpha*(1+Forcing(j+1))*(1/(1+incz^n))-incx;
    y2=incx-incy;
    z2=incy-incz;

    x(j+1)=x(j)+(x1+x2)*time_step/2;
    y(j+1)=y(j)+(y1+y2)*time_step/2;
    z(j+1)=z(j)+(z1+z2)*time_step/2;
end
% T_shift=60;


%% spike calculation

[val1,locs1] = findpeaks(x,'MinPeakProminence',(max(x)-min(x))/20);
[val2,locs2] = findpeaks(Forcing);
%  plot(time(locs1),val1,'bs');
%  hold on
%  plot(time(locs2),val2,'ks');
period_system=time(locs1(end))-time(locs1(end-1));
period_forcing=time(locs2(end))-time(locs2(end-1));
period_ratio=period_system/period_forcing;
disp(['spike period: ', num2str(period_system)]);
disp(['Forcing period: ', num2str(period_forcing)]);
disp(['Period zRatio: ', num2str(period_ratio)]);

% num_spikes = numel(locs1);
% start_idx = floor(2*length(time)/3);
start_idx = T_shift*100;
locs1_in_window = locs1(locs1 >= start_idx);  % only peaks in the last 1/3
num_spikes = numel(locs1_in_window);

% num_pulses = floor(Total_time / per);
num_pulses = floor((Total_time - T_shift) / per);

ratio=num_spikes/num_pulses;
disp(['Number of spikes within the interval: ', num2str(num_spikes)]);
disp(['Number of pulses within the interval: ', num2str(num_pulses)]);
disp(['Ratio: ', num2str(ratio)]);


%%
figure(1)
subplot(4,1,[1 2 3])
plot(time-T_shift,x,'b','linewidth',3);

xticklabels([]); % Remove x-axis labels in this subplot
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18);
ylabel('$X$', 'Interpreter', 'latex')
ytickangle(0)
xlim([0 Total_time-T_shift])
ylim_max=max(x);
ylim_min=min(x);
ylim([ylim_min-0.1 ylim_max+0.1]);
% ylim([1.2 2.3]);

subplot(4,1,4);
plot(time-T_shift,Forcing,'r','linewidth',3);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18);
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$F$', 'Interpreter', 'latex')
xlim([0 Total_time-T_shift])
ylim([0 Amp+0.11])


