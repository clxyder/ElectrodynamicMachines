%% Project 2 - Arbitrary Reference Frame
% Carlos A. Wong
% Stephanie Damas

%% Overview
% Three-phase quantities are transformed to dqz quantities for further
% analysis.
%% Inputs:
% 
%% Outputs:
% 
%% Setup
clear
close all;
fntSz = 11;
%% Parameters
v_s = 10/sqrt(2);       % V         source
r_s = 0.216;            % ohms      series resistance
w_e = 377;              % rad/s     angular velocity
wL_es = 1.09;           % ohms
L_s = wL_es/w_e;        % H         inductance

z_s = r_s + 1j*wL_es;

tau = L_s/r_s;
alpha = atan(wL_es/r_s);

ramp_start = 0.05;

test_case = 1;
scale = 1;

%% Simulation
simTime = 0.12;
% sim('proj2.slx',simTime)

%% Unbalancing of System
system = {{'balanced' 1}; {'unbalanced' 0.6}};

%% Figures
plotRow = 6;
plotRow2 = 3;
plotCol = 1;
xmax = 1200;
ymax = 12;
pos = [4.2 1.2 7.17 7];
labels = {'v_{qs} (V)'; 'i_{qs} (A)'; 
    'v_{ds} (V)'; 'i_{ds} (A)'; 
    'Power (W)'; '\omega (rad/s)'};
labels2 = {'v_{as} (V)'; 'v_{bs} (V)'; 'v_{cs} (V)';
           'i_{as} (A)'; 'i_{bs} (A)'; 'i_{cs} (A)';};
names = {'Figure_3_9_1' 'Figure_3_9_2' 'Figure_3_9_3'};
modes = {'Stationary' 'Synchronous', 'Jump & Run'};
letters = lower(char((1:numel(labels)) + 64));

%% All Scenarios
for kk = 1:length(system)
    config = system{kk}{1};
    scale =  system{kk}{2};
    
    %% Source Figures
    test_case = 1;sim('proj2.slx',simTime)
    data2 = {v_abc.Data(:,1); v_abc.Data(:,2); v_abc.Data(:,3);
         i_abc.Data(:,1); i_abc.Data(:,2); i_abc.Data(:,3)};
    H2 = figure(kk+10);
    H2.Name = [config '_source'];
    for ii = 1:plotRow
        subplot(plotRow, plotCol, ii)
        plot(data2{ii}); xlim([0 xmax]); grid on; box on
        ylabel(labels2{ii}, 'FontWeight','bold', ...
            'interpreter','tex','FontSize', fntSz)
        set(gca,'FontSize',fntSz)
        letter = ['(' letters(ii) ')'];
        annotation('textbox', [0.91, 0.93-0.1422*(ii-1), 0, 0], ...
            'string', letter,'FontWeight', 'bold','FontSize',fntSz)
        if ii == 6
            xlabel('Time (ms)','FontWeight','bold','FontSize',fntSz)
        end
    end
    
    set(H2, 'Units', 'in', 'Position', pos, ...
        'PaperPositionMode','auto','PaperSize', 0.9*pos(3:4))
    saveas(H2,['figures/' H2.Name], 'pdf')

    for jj = 1:length(names)
        test_case = jj;
        sim('proj2.slx',simTime)
        data = {v_qdz.Data(:,1); i_qdz.Data(:,1);
                v_qdz.Data(:,2); i_qdz.Data(:,2);
                p_qdz.Data(:,1); omega.Data(:,1)};

        name = names{jj};
        %% Book Figures
        if kk == 1
            H = figure(jj);
        else
            H = figure(jj+4);
        end
        H.Name = [config '_' name];
        for ii = 1:plotRow
            subplot(plotRow, plotCol, ii)
            plot(data{ii}); xlim([0 xmax]); grid on; box on
            ylabel(labels{ii}, 'FontWeight','bold', ...
                'interpreter','tex','FontSize', fntSz)
            set(gca,'FontSize',fntSz)
            letter = ['(' letters(ii) ')'];
            annotation('textbox', [0.91, 0.93-0.1422*(ii-1), 0, 0], ...
                'string', letter,'FontWeight', 'bold','FontSize',fntSz)
            if ii == 6
                xlabel('Time (ms)','FontWeight','bold','FontSize',fntSz)
            end
            if ii == 1 || ii == 3
                ylim([-ymax ymax])
            end
        end
        set(H, 'Units', 'in', 'Position', pos, ...
            'PaperPositionMode','auto','PaperSize', 0.9*pos(3:4))
        saveas(H,['figures/' H.Name], 'pdf')
        
        %% Three Phase Power Figures
        h = figure(kk*2);
        h.Name = [config '_abc_power'];
        subplot(plotRow2,plotCol, jj); plot(p_abc.Data(:,1))
        ylabel('Power (W)','FontWeight','bold','FontSize', fntSz)
        xlim([0 xmax]); grid on; box on; set(gca,'FontSize',fntSz)
        text(900, 90, modes{jj},'FontWeight','bold','FontSize', fntSz)
        letter = ['(' letters(jj) ')'];
        annotation('textbox', [0.91, 0.93-0.3*(jj-1), 0, 0], ...
                'string',letter,'FontWeight','bold','FontSize',fntSz)
        if jj == 3
            xlabel('Time (ms)','FontWeight','bold','FontSize',fntSz)
        end
        h2 = figure(kk+15);
        h2.Name = [config '_abc_current'];
        subplot(plotRow2,plotCol, jj); hold on;
        plot(i_abc.Data(:,1)); plot(i_abc_inv.Data(:,1)-1.5,'-.');hold off;
        legend('Before transformation', 'After inverse transformation')
        ylabel('Current (A)','FontWeight','bold','FontSize', fntSz)
        xlim([0 xmax]); grid on; box on; set(gca,'FontSize',fntSz)
        text(900, 90, modes{jj},'FontWeight','bold','FontSize', fntSz)
        letter = ['(' letters(jj) ')'];
        annotation('textbox', [0.91, 0.93-0.3*(jj-1), 0, 0], ...
                'string',letter,'FontWeight','bold','FontSize',fntSz)
        if jj == 3
            xlabel('Time (ms)','FontWeight','bold','FontSize',fntSz)
        end
    end
    pos2 = h.PaperPosition;
    set(h, 'PaperSize', 0.97*pos2(3:4))
    saveas(gcf,['figures/' h.Name], 'pdf')
    pos2 = h2.PaperPosition;
    set(h2, 'PaperSize', 0.97*pos2(3:4))
    saveas(h2,['figures/' h2.Name], 'pdf')
end