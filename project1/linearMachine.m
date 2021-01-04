%% Function Description
% Given the resistance, inductance and mass plot the:
%       e_f         electromotive force
%       lambda      flux linkage
%       i           current
%       f_e         electromagnetic force
%       x           position
%       W_e         electrical energy
%       W_f         coupling field energy
%       W_m         mechanical energy
%% Inputs
%       R           resistance
%       L           inductance
%       M           mass
%% Outputs
%       v           voltage
%% Static Parameters
R   = 10;         % ohms
x0  = 3e-3;       % millimeters         (Initial Position)
K   = 2667;       % N/m                 (Spring Coefficeint)
k   = 6.293e-5;   % H*m                 (Inductance per meter)
v_o = 5;          % V                   (Voltage)
F_o = 4;          % N                   (Force)
% Dynamic Parameters
l   = 0;          % H                   (Line Inductance)
M   = 0.055;      % kg                  (Mass of Moving Member)
D   = 4;          % N*s/m               (Damping Coefficient)
%% Pulsed Waveforms
% v_o pulse T = 0.9, duty cycle = 90, delay = .12
% F_o pulse T = 0.75, duty cycle = 60, delay = 0.27
%% Varying Parameters
scale =  0.1;
% R = scale*R;            % Constants: F_o, v_o, M
% F_o = scale*F_o;        % Constants: R, v_o, M
% M = scale*M;            % Constants: F_o, v_o, R

%% SimuLink Model
Simulation_Time = 1;
sim('LinearMachineSim',Simulation_Time)

%% Plot Generation
figure = {'Figure 1.3-10'
        'Figure 1.3-11'
        ['Varying Resistance R = ' num2str(R) ' ohm']
        ['Varying Force F = ' num2str(F_o) ' N']
        ['Varying Mass M = ' num2str(M) ' kg']};
plotGraph = 2;

if plotGraph == 1
    plotRow = 4; plotCol = 2;
    figure; sgtitle(figure(2))
    %% Electromotive Force
    subplot(plotRow,plotCol,1); plot(emf)
    title('Electromotive Force')
    ylabel('Voltage (V)'); xlabel('Time (s)');
    %% Flux Linkage
    subplot(plotRow,plotCol,2); plot(flux_linkage)
    title('Flux Linkage')
    ylabel('Volt-Seconds (V.s)'); xlabel('Time (s)');
    %% Current
    subplot(plotRow,plotCol,3); plot(current)
    title('Current')
    ylabel('Current (A)'); xlabel('Time (s)');
    %% Electromagnetic Force
    subplot(plotRow,plotCol,4); plot(f_e)
    title('Electromagnetic Force')
    ylabel('Newton (N)'); xlabel('Time (s)');
    %% Position
    subplot(plotRow,plotCol,5); plot(position)
    title('Position')
    ylabel('Distance (mm)'); xlabel('Time (s)');
    %% Electrical Energy
    subplot(plotRow,plotCol,6); plot(elec_energy)
    title('Electrical Energy')
    ylabel('Joules (mJ)'); xlabel('Time (s)');
    %% Field Energy
    subplot(plotRow,plotCol,7); plot(field_energy)
    title('Field Energy')
    ylabel('Joules (mJ)'); xlabel('Time (s)');
    %% Mechanical Energy
    subplot(plotRow,plotCol,8); plot(mech_energy)
    title('Mechanical Energy')
    ylabel('Joules (mJ)'); xlabel('Time (s)');
elseif plotGraph == 2
    plotRow = 2; plotCol = 2;
    figure; sgtitle(figure(5))
    %% Electromotive Force
    subplot(plotRow,plotCol,1); plot(emf)
    title('Electromotive Force')
    ylabel('Voltage (V)'); xlabel('Time (s)');
    %% Current
    subplot(plotRow,plotCol,2); plot(current)
    title('Current')
    ylabel('Current (A)'); xlabel('Time (s)');
    %% Electromagnetic Force
    subplot(plotRow,plotCol,3); plot(f_e)
    title('Electromagnetic Force')
    ylabel('Newton (N)'); xlabel('Time (s)');
    %% Velocity
    subplot(plotRow,plotCol,4); plot(velocity)
    title('Velocity')
    ylabel('Speed (mm/s)'); xlabel('Time (s)');
end