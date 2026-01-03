clear all;
close all;
clc

%% Fixed Parameters
% vm=1; 
km=0.1; vp=0.5;
kp1=10; kp2=0.03; kp3=0.1;
Pcrit=0.1; Jp=0.05;
Keq=200;

vm=0.2207; 
A=0.01;
% Keq=3.3;

%% Time setup
tstart=0;
tIni=0;
dt=0.01;

%% Arrays of PHP, T_left, T_right (same length)
php_array = [0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.92 0.922	0.93 0.932];
T_left_array = [58.5	58	57.9	57.85	58	58 57.75	57.5	57.75	57.75	58.3	58.85];
T_right_array = [58.5	58.8	60.2	60.75	60.75 61.2	61	60.25	59.75	59.65	59.25	58.85];
% % % php_array = [0.01, 0.1];
% % % T_left_array = [24.8, 25];
% % % T_right_array = [26.3, 26];
%% Initialize master results array
all_results = [];

% Initialize cell array
results_cell = cell(length(php_array), 1); % because it is not recommended to do inside the loop in parallel

%% Loop over parameter sets
parfor idx = 1:length(php_array)
    fprintf("Running idx = %d on worker %d\n", idx, getCurrentTask().ID); % to show how parallel work

    php = php_array(idx);
    T_left = T_left_array(idx);
    T_right = T_right_array(idx);

    % Create 21 per values
    per_values = linspace(T_left, T_right, 21);

    % Preallocate
    fourth_peak_times = NaN(length(per_values), 1);
    fourth_pulse_start_times = NaN(length(per_values), 1);
    phase = NaN(length(per_values), 1);

    for i = 1:length(per_values)
        per = per_values(i);
        tend = per * 10;
        tspan = tIni:dt:tend;

        IC = [0.1 0.1];
        options=odeset('AbsTol',1e-5,'relTol',1e-5,'InitialStep',1e-3,'MaxStep',0.01);
        [t1,y1]=ode15s(@NTdy,tspan,IC,options, ...
            vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per);
        IC = [y1(end,1) y1(end,2)];
        [ti1,yy1]=ode15s(@NTdy,tspan,IC,options, ...
            vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per);

        M = yy1(:,1);

        % Counting peaks
        [pks, locs] = findpeaks(M,'MinPeakProminence',(max(M)-min(M))/10);
        valid_locs = locs(locs >= tstart*100);

        %% Extract the time and the value of the 4th peak of M
        if numel(valid_locs) >= 4
            fourth_peak_times(i) = valid_locs(4) * dt;
        end

        %% Find the starting point of the 3rd pulse of F
        fourth_pulse_start_times(i) = 3 * per;

        %% Calculate the phase (time_{oscillator} - time_{forcing})/period
        if ~isnan(fourth_peak_times(i))
            phase(i) = mod((fourth_peak_times(i) - fourth_pulse_start_times(i)) / per, 1);
        end
    end

    % Store current results with php as the first column
    block_results = [repmat(php, length(per_values), 1), ...
        per_values(:), ...
        fourth_peak_times(:), ...
        fourth_pulse_start_times(:), ...
        phase(:)];

    %     all_results = [all_results; block_results];  % Append to the full result
    results_cell{idx} = block_results;

end
all_results = vertcat(results_cell{:});  % Combine safely after parallel run


%% Save to a single CSV file
% file_name = sprintf('Keq_%.2f_A_%.2f_phase_results3.csv', Keq, A);
file_name = sprintf('Vm_%.2f_A_%.2f_phase_results4.csv', vm, A);

csvwrite(file_name, all_results);

%% ODE function
function dy = NTdy(ti1,y,vm, km, vp,kp1,kp2,kp3,Keq,Pcrit,Jp,A,php,per)
F = A*(1 - heaviside(mod(ti1,per) - per*(1-php)));
q = 2 / (1 + sqrt(1 + 8*Keq*y(2)));
dy(1,1) = (1 + F)*vm / (1 + (y(2)*(1 - q)/(2*Pcrit))^2) - km*y(1);
dy(2,1) = vp*y(1) - (kp1*y(2)*q + kp2*y(2)) / (Jp + y(2)) - kp3*y(2);
end
