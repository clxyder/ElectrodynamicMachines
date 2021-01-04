%% Project 4 - dq0 model & steady – state model
% Carlos A. Wong & Taylor T. Stamm
close all
%% Electrical Parameters for 2.23 kW - 3 HP
Name        = '3_HP';   % Type of Machine
v_rms       = 220;              % v rms line to line (V)
v_ll        = v_rms/sqrt(3);    % Peak Line to Line Voltage (V)
v_pk        = v_rms*sqrt(2/3);  % Peak Line to Ground Voltage (V)
f           = 60;               % frequency (Hz)
P           = 4;        % number of poles (#)
wb          = 2*pi*f;   % base speed (rad/s)
we          = wb;       % excitation frequency (rad/s)
B           = 0.00001;  % Damping Coefficient
rs          = 0.45;     % Stator Winding Resistance (ohms) 
X_lr_p      = 0.75;     % Referred Rotor Leakage Reactance (ohms)
rrp         = 0.8;      % Referred Rotor Winding Resistance (ohms)
X_ls        = 0.75;     % Stator Leakage Reactance (ohms)
X_M         = 27;       % Magnetized Mutual Reactance (ohms)
J           = 0.09;     % Inertia (kg-m^2)
w_ref       = 0;        % reference frame speed (rad/s)
to_rpm      = 30/pi;    % rad/s to rpm conversion
sim_time    = 5;        % Simulation Time (s)
load_time   = 0.5;
% load_time   = 6;

%% Simmulation Params
params = [P wb X_lr_p X_ls X_M rs rrp B J];

%% Simulating
sim('steady_state',sim_time);

%% Data Extraction
t = Te.Time;
torque = Te.Data;
speed = wr.Data;
w_rpm = to_rpm*speed;
va = v_abc.Data(:,1);
vb = v_abc.Data(:,2);
vc = v_abc.Data(:,3);
ias = i_abcs.Data(:,1);
ibs = i_abcs.Data(:,2);
ics = i_abcs.Data(:,3);
iarp = i_abcrp.Data(:,1);
ibrp = i_abcrp.Data(:,2);
icrp = i_abcrp.Data(:,3);
p_in = (va.*ias+vb.*ibs+vc.*ics);
p_shaft = (2/P)*speed.*torque;
%% Stall Torque
k = find(speed<1e-6,3,'first');
t_stall = torque(k);
disp(['Stall Torque = ' num2str(t_stall)])
disp(['Stall Torque @ 25% = ' num2str(0.25*t_stall)])

%% Torque
H(1) = figure('Name','Torque');
plot(t,torque,'black', 'LineWidth',1.5)
xlabel('Time (s)')
ylabel('\tau_e (N-m)')
set(gca,'FontSize',16);
title('Te vs. Time')
legend('Te')
H(1).Units = 'in';
H(1).PaperSize = 0.99*H(1).Position(3:4);
saveas(H(1),'ssTeTime','pdf')
%% Speed
% figure('Name','Speed');
% plot(t,speed)
%% Speed vs Torque
H(2) = figure('Name','Speed vs Torque');
plot(w_rpm,torque,'black', 'LineWidth',1.5)
xlabel('\omega_r (rpm)')
ylabel('\tau_e (N-m)')
set(gca,'FontSize',16);
title('Te Vs. Rotor Speed')
legend('Te')
H(2).Units = 'in';
H(2).PaperSize = 0.99*H(2).Position(3:4);
saveas(H(2),'ssSpeedTorque','pdf')
%% Input Power
H(3) = figure('Name','Input Power');
plot(t,p_in,'black', 'LineWidth',1.5)
xlim([0.5 5])
xlabel('Time (s)')
ylabel('Power (W)')
legend('P_{input}')
title('Input Power vs. Time')
set(gca,'FontSize',16);
H(3).Units = 'in';
H(3).PaperSize = 0.99*H(3).Position(3:4);
saveas(H(3),'ssPowerIn','pdf')
%% Shaft Power
H(4) = figure('Name','Shaft Power');
plot(t,p_shaft,'black', 'LineWidth',1.5)
xlim([0.5 5])
xlabel('Time (s)')
ylabel('Power (W)')
legend('P_{shaft}')
title('Shaft Power vs. Time')
set(gca,'FontSize',16);
H(4).Units = 'in';
H(4).PaperSize = 0.99*H(4).Position(3:4);
saveas(H(4),'ssPowerShaft','pdf')