%% Project 3 - Induction Machine Simulation
% Carlos A. Wong
% Hebert Lopez
%% Inputs
% v_rms         RMS Voltage Line to Line (V)
% P             Number of Magnetic Poles (#)
% w_b           Base Speed (rad/s)
% r_s           Stator Winding Resistance (ohms)
% X_lr_p        Referred Rotor Leakage Reactance (ohms)
% r_rp          Referred Rotor Winding Resistance (ohms)
% X_ls          Stator Leakage Reactance (ohms)
% X_M           Magnetized Mutual Reactance (ohms)
% J             Inertia (kg-m^2)
%% Outputs

%% Figure Configurations
clear;close all;clc;
warning('off','MATLAB:print:FigureTooLargeForPage')
fig = 0;fig_save = 0;
fig_sel = questdlg('Figure option?', ...
    'Process Figures','Figures and Save','Figures','None','None');
if strcmp(fig_sel, 'Figures and Save')
    fig = 1; fig_save = 1;
elseif strcmp(fig_sel, 'Figures')
    fig = 1; fig_save = 0;
end
fntSz = 11;
labels = {'i_{as} (A)';'i_{bs} (A)';'i_{cs} (A)';
    'i_{ar}'' (A)';'i_{br}'' (A)';'i_{cr}'' (A)';
    '\tau_{e} (N-m)'};
letters = lower(char((1:numel(labels)) + 64));
%% Electrical Parameters
v_rms   = 220;              % v rms line to line (V)
v_pk    = v_rms*sqrt(2/3);  % Peak Line to Ground Voltage (V)
P       = 4;                % number of poles (#)
f       = 60;               % frequency (Hz)
w_b     = 2*pi*f;           % base speed (rad/s)
B       = 0.00001;          % Damping Coefficient
cnv     = 0.7457;           % kW = cnv * HP

%% 2.23 kW - 3 HP
hp3 = struct;
hp3.Name        = '3_HP';   % Type of Machine
hp3.Rating      = 2.23e3;   % Power Rating (W)
hp3.r_s         = 0.45;     % Stator Winding Resistance (ohms)
hp3.X_lr_p      = 0.75;     % Referred Rotor Leakage Reactance (ohms)
hp3.r_rp        = 0.8;      % Referred Rotor Winding Resistance (ohms)
hp3.X_ls        = 0.75;     % Stator Leakage Reactance (ohms)
hp3.X_M         = 27;       % Magnetized Mutual Reactance (ohms)
hp3.J           = 0.09;     % Inertia (kg-m^2)
hp3.sim_time    = 0.6;      % Simulation Time (s)
hp3.tol         = 0.3;      % Steady State Tolerance ()
hp3.load_time   = 0.4;      % Time to Load machine (s)

%% 5.22 kW - 7 HP
hp7 = struct;
hp7.Name        = '7_HP';   % Type of Machine
hp7.Rating      = 5.22e3;   % Power Rating (W)
hp7.r_s         = 0.3;      % Stator Winding Resistance (ohms)
hp7.X_lr_p      = 0.27;     % Referred Rotor Leakage Reactance (ohms)
hp7.r_rp        = 0.15;     % Referred Rotor Winding Resistance (ohms)
hp7.X_ls        = 0.57;     % Stator Leakage Reactance (ohms)
hp7.X_M         = 20;       % Magnetized Mutual Reactance (ohms)
hp7.J           = 0.25;     % Inertia (kg-m^2)
hp7.sim_time    = 0.9;      % Simulation Time (s)
hp7.tol         = 5e-3;     % Steady State Tolerance ()
hp7.load_time   = 0.72;     % Time to Load machine (s)

%% Test Selection
dut_sel = questdlg('Machine to simulate?', ...
    'Test Selector',hp3.Name,hp7.Name,'Both',hp3.Name);

if strcmp(dut_sel,hp3.Name)
    duts = {hp3};
elseif strcmp(dut_sel,hp7.Name)
    duts = {hp7};
elseif strcmp(dut_sel,'Both')
    duts = {hp3 hp7};
end

%% Iterate Test Simulations
for ii = 1:length(duts)
dut = duts{ii};
%% Calculation of Machine Specific Paramters
L_lr_p      = dut.X_lr_p/w_b;       % Referred Rotor Leakage Inductance (L)
L_ls        = dut.X_ls/w_b;         % Leakage Inductance from Stator (H)
L_ms        = dut.X_M/w_b;          % Magnetized Inductance(H)
r_s         = dut.r_s;
r_rp        = dut.r_rp;
J           = dut.J;
params = [L_ls L_ms L_lr_p P r_s r_rp J B];

%% Modify Simulation
vs_toff = 3;        % time to turn of stator voltage supply (s)
load = 0;           % load torque (N-m)
if load
load_tL = dut.load_time;        % time to load machine (s)
else
load_tL = 3;
end
%% Simulation
disp('-----------------------')
disp(['Device under Test: ' dut.Name])
disp('-----------------------')
sim('induction_machine',dut.sim_time);

%% Prelimanary Analysis
% Stator and Rotor Currents with Torque vs. Time
if fig
H(1) = figure('Name',[dut.Name '_Currents_Torque'],'NumberTitle','off');
H(1).Position = [0.15 0.6 2.4 2.7]*300;
plotData = {i_s.Data(:,1); i_s.Data(:,2); i_s.Data(:,3);
        i_rp.Data(:,1);i_rp.Data(:,2);i_rp.Data(:,3);
        t_e.Data};
for jj = 1:length(plotData)
    subplot(7,1,jj);plot(tout,plotData{jj});grid on; box on;
    ylabel(labels{jj},'fontweight','bold','FontSize',fntSz);
    annotation('textbox', [0.91, 0.93-0.1215*(jj-1), 0, 0], ...
        'string', ['(' letters(jj) ')'],'FontWeight', 'bold','FontSize',fntSz)
    set(gca,'FontSize',fntSz+2.1)
end
xlabel('Time (seconds)','fontweight','bold','FontSize',fntSz)

H(1).Units = 'in';
H(1).PaperSize = 0.96*H(1).Position(3:4);
if fig_save
    saveas(H(1),['figures/' H(1).Name],'pdf')
end
% Torque vs.Speed
H(2) = figure('Name',[dut.Name '_Torque_vs_Speed'],'NumberTitle','off');
H(2).Position = [690 180 H(2).Position(3:4)];
plot(w_r.Data, t_e.Data);grid on; box on;
ylabel('\tau_{e} (N-m)','fontweight','bold','FontSize',fntSz)
xlabel('\omega_{r} (rad/s)','fontweight','bold','FontSize',fntSz)
set(gca,'FontSize',fntSz)

H(2).Units = 'in';
H(2).PaperSize = 0.99*H(2).Position(3:4);
if fig_save
    saveas(H(2),['figures/' H(2).Name],'pdf')
end
%% Advanced Analysis
H(3) = figure('Name',[dut.Name '_Power'],'NumberTitle','off');
H(3).Position = [1260 180 H(3).Position(3:4)];
plot(tout,[p_stator.Data,p_rotor.Data,p_shaft.Data])
ylabel('Power (W)','fontweight','bold','FontSize',fntSz)
xlabel('Time (seconds)','fontweight','bold','FontSize',fntSz)
grid on; box on;
legend('P_{stator}','P_{rotor}','P_{shaft}')
ax = gca;
ax.YAxis.Exponent = 3;
set(gca,'FontSize',fntSz)

H(3).Units = 'in';
H(3).PaperSize = 0.99*H(3).Position(3:4);
if fig_save
    saveas(H(3),['figures/' H(3).Name],'pdf')
end
end
%% Peak Inrush Current
i_pkinrush = max(max(abs(i_s.Data)));
disp(['Peak Inrush Current: ' num2str(i_pkinrush) ' (A)'])
%% Peak to Peak Torque
te_pkpk = [min(t_e.Data) max(t_e.Data)]; % [min_pk max_pk]
disp('Peak to Peak Torque:')
disp(['min: ' num2str(te_pkpk(1)) '(N-m)'])
disp(['max: ' num2str(te_pkpk(2)) ' (N-m)'])
%% Time to reach Steady State
dwr = gradient(w_r.Data);
steady = find(abs(dwr)<dut.tol,1,'last');
% steady = find(abs(dw_r.Data)<1,1,'last');
% figure;plot(tout,[dwr,zeros(size(tout))])
% figure(9+ii);plot(tout,[dw_r.Data,zeros(size(tout))])
t_steady = tout(steady);
disp(['Time to reach Steady State: ' num2str(t_steady) ' (s)'])
%% Overshoot of Steady State
steady_te = t_e.Data(steady);
disp(['Torque at Steady State: ' num2str(steady_te) ' (N-m)'])
ovr_shoot = 100*abs(te_pkpk(2)-steady_te)/steady_te;
disp(['Percent Overshoot at Steady State: ' num2str(ovr_shoot) ' (%)'])
end