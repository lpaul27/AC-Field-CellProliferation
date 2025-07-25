function [] = Visualize(x_time,y_time, theta_time, time_control, dispAvg, directedness, fields, cell_posData)
%% Function for static visualization
% Optimal for large sets of data

%% Plots:
% Boolean logic for on and off
%% - Displacement:
% Tracks cell position in a 2D x-y graph over time
% **displacement = 1**: ON
% **displacement = 0**: OFF
%% - (2D) Directionality:
% Tracks cell direction of propogation over time
% **directionality = 1**: ON
% **directionality = 0**: OFF
%% - Polar histogram:
% Tracks distribution of cell direction across a polar domain
% **polarhist = 1**: ON
% **polarhist = 0**: OFF
%% - (1D) Directionality:
% Directionality between a single axis of movement (tiled)
% **dim1directionality = 1**: ON
% **dim1directionality = 0**: OFF

%% Begin Function
global NumCells runTime displacement live dim2directionality polarhist dim1directionality dim1displacement...         %#ok<GVMIS>
    disphist directednessplot ExMax displacement3by2 runs densityplotDIR densityplotDISP %#ok<GVMIS>

fieldboundslx = [-150, -200, -200, -250, -250, -250];
fieldboundshx = [150, 100, 100, 50, 50, 50];
fieldboundshy = [150, 150, 150,  150, 150, 150];
fieldboundsly = [-150, -150, -150, -150, -150, -150];
% runs nothing if live simulation is on

if(~live)
    %% Displacement track graph
    % Visualization of cell proliferation

    if(displacement)

        %calculations
        x_avg = mean((x_time - x_time(1,:)), 2);
        y_avg = mean((y_time - y_time(1,:)), 2);

        %plotting
        figure
        hold on;
        plot((x_time - x_time(1,:)), (y_time - y_time(1,:)))
        hold on;
        plot((x_avg - x_avg(1,1)), (y_avg - y_avg(1,1)), 'k', 'LineWidth', 3);
        xline(0, '-');
        yline(0, '-');
        xlabel('x-displacement (a.u)')
        ylabel('y-displacement (a.u)')
        print1 = sprintf('%.0f mV/mm', (ExMax * 3750));
        lgd = legend(print1);
        title(lgd,'EF strength');
        xlim([-200, 200]);
        ylim([-200,200]);
    end

    %% 1D Directionality Graph
    % Visualization of allignment to 180 degree field change

    quantitative_data = 1;

    if(dim1directionality)
        cosTheta = zeros(runTime, runs);
        sumCellAnglex = zeros(runTime, runs);
        directionalityX = zeros(runTime, runs);

        for i = 1:runs 
            cosTheta(:,i) = mean((cell_posData(i).direct5),2);
            sumCellAnglex(:,i) = sum((cosTheta(:,i)),2);
            directionalityX(:,i) = sumCellAnglex(:,i);
        end
        directionalityMean4 = mean(cosTheta,2);

        if(quantitative_data)
            confIntXrise = find(directionalityX >= (max(directionalityX) - 0.05));
            confIntXfall = find(directionalityX <= (min(directionalityX(10:end)) + 0.05) & gradient(directionalityX) < 0);
            XriseTime = confIntXrise(1,1);
            XfallTime = confIntXfall(1,1);
        end
        figure
        for i = 1:runs
            scatter(2*time_control, cosTheta(:,i));
            hold on
        end
        hold on
        plot((time_control * 2), directionalityMean4, 'Linewidth', 2, 'Color', [0 0 0])
        xlabel('Time (minutes)', 'FontSize', 18);  ylabel('Directionality (\Phi_{x})', 'FontSize', 18);
        %xticks = 2*(1:runTime)';
        xticks([0 60 120 180 240])
        xticklabels([0 60 120 180 240]);
        xline((runTime / 2) * 2,'-.', 'FLIP', 'LineWidth',2, 'FontSize', 14);
        ylim([-1.0,1.0]); xlim([0, runTime * 2]);
        
        if(quantitative_data)
            print1 = sprintf('X Rise: %f', XriseTime);
            print2 = sprintf('X Fall: %f', (XfallTime - (runTime / 2)));
            xline(XriseTime, ':', print1, 'Color', 'b', 'LineWidth',1);
            xline(XfallTime, ':', print2, 'Color', 'b', 'LineWidth',1);
        end

        yline(0);
        yline(1, '--');
        yline(-1, '--');
        txt1 = {'FIELD RIGHT'}; txt2 = {'FIELD LEFT'};
        text((runTime*2 /9), 1.1, txt1, 'FontSize',14);
        text((5*runTime /4), 1.1, txt2, 'FontSize',14);
        set(gcf,'Color','white');
        set(gca, 'FontSize', 18);
        exportgraphics(gcf, 'direcRise.pdf', 'ContentType', 'vector', 'Resolution',1000);    
    end
    %% 2D Directionality Graph
    % Visualization of allignent to a direction

    if(dim2directionality)

        %calculations
        cosTheta = cos(theta_time);
        sinTheta = sin(theta_time);
        sumCellAnglex = sum(cosTheta,2);
        sumCellAngley = sum(sinTheta,2);
        directionalityX = sumCellAnglex ./ NumCells;
        directionalityY = sumCellAngley ./ NumCells;
        confIntXrise = find(directionalityX >= (max(directionalityX) - 0.05));
        confIntXfall = find(directionalityX <= (min(directionalityX(10:end)) + 0.05) & gradient(directionalityX) < 0);
        confIntYrise = find(directionalityY >= (max(directionalityY) - 0.05));
        % lead point is the first point found
        XriseTime = confIntXrise(1,1);
        YriseTime = confIntYrise(1,1) ;
        XfallTime = confIntXfall(2,1);
        % Plotting
        figure
        plot(time_control, directionalityX)
        hold on;
        plot(time_control, directionalityY, '--');
        xlabel('Time (steps)');  ylabel('Directionality (\Phi)');
        ylim([-0.2,1.2]); xlim([0, runTime]);
        xline((runTime / 2),'-.', 'TURN', 'LineWidth',2);
        print1 = sprintf('X Rise: %f', XriseTime);
        print2 = sprintf('X Fall: %f', (XfallTime - (runTime / 2)));
        print3 = sprintf('Y Rise: %f', (YriseTime - (runTime / 2)));
        xline(XriseTime, ':', print1, 'Color', 'b', 'LineWidth',1);
        xline(XfallTime, ':', print2, 'Color', 'b', 'LineWidth',1);
        xline(YriseTime, ':', print3, 'Color','r', 'LineWidth',1);

        legend('$\mathrm{\Phi_{x}}$', '$\mathrm{\Phi_{y}}$','Interpreter', 'latex','Location', 'southeast');
        txt1 = {'FIELD RIGHT'}; txt2 = {'FIELD UP'};
        text((runTime /8), 1.1, txt1);
        text((5*runTime /8), 1.1, txt2);
    end
    %% Distribution Plot visualization
    % Visualization of directional distribution

    if(polarhist)

        % NOT FINISHED!
        % !!!!!!!!!!!!! %

        %calculations
        stat_raw = reshape(theta_time, 1, []);
        stat_trajAvg = mean(theta_time, "all");

        %plotting
        polarhistogram(stat_raw', 36);
        hold on;
        rlim_vals = rlim;
        vector_length = rlim_vals(2);
        polarplot([stat_trajAvg stat_trajAvg], [0, vector_length], 'r-', 'LineWidth',2);
    end
    %% 1D Directionality time based plots
    % For visualization of directionality in time represented across a 1D axis

    if(dim1displacement)
        tiledlayout(2,1)
        nexttile
        plot((x_time(1:floor(runTime/2), :) - x_time(1,:)), (y_time(1: floor(runTime / 2), :) - y_time(1,:)))
        title('Field Direction: +x')

        nexttile
        plot((x_time(floor(runTime/2) + 1: runTime, :) - x_time(runTime / 2,:)), (y_time(floor(runTime / 2) + 1: runTime, :) - y_time(runTime/ 2 + 1,:)))
        title('Field Direction: -x')
    end % end single axis directionality

    %% distribution of data
    if(disphist)

        groups = {'No EF', '30mV/mm', '50mV/mm', '75mV/mm', '100mV/mm', '200mV/mm'};

        % Example data
        data_sim = dispAvg(:, 1);  % Make sure dispAvg is a column vector
        data_exp = [0.45; 0.8; 1.35; 1.3; 1.45; 1.6];
        electric_fields = [0,30, 50, 75, 100, 200];

        % Combine into matrix (columns = datasets, rows = groups)
        data = [data_sim, data_exp];  % [6x2] matrix

        % Create grouped bar chart
        figure;
        b = bar(electric_fields, data, 'grouped');  % grouped bar by default

        % Set colors
        b(1).FaceColor = [0 0.5 1];   % blue
        b(2).FaceColor = [1 0.4 0];   % orange

        % Label settings
        set(gca, 'XTickLabel', electric_fields, 'XTickLabel', groups, 'XTickLabelRotation', 45);
        ylabel('Displacement speed (\mum/min)');
        ylim([0 3]);

        legend({'Simulation', 'Experiment'}, 'Location', 'northwest');
        box on;
        set(gca,'fontsize',14);
    end
    if(directednessplot)
        groups = {'No EF', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
        data_tmp = zeros(runTime,1);
        err_tmp = zeros(runTime,1);
        % Example data
        for i = 1:length(groups)
            
            data_tmp(:, i) = cell_posData(i).directionality_mean;  % Make sure dispAvg is a column vector
            err_tmp(:,i) = cell_posData(i).directionalitySD;
            data_sim = data_tmp(runTime,:)';
            err = err_tmp(runTime,:)';
        end
        if(quantitative_data)
            data_exp = [0.05; 0.4; 0.5; 0.6; 0.8; 0.9];
            electric_fields = [0,30, 50, 75, 100, 200];
        end
        % Combine into matrix (columns = datasets, rows = groups)
        if(quantitative_data)
            data = [data_sim, data_exp];  % [6x2] matrix
        else
            data = data_sim;
        end
        % Create grouped bar chart
        figure;
        b = bar(data);
        hold on
        er = errorbar((1:length(groups)),data, err, 'o');
        er.Color = [0 0 0];

        % Set colors (optional)
        b(1).FaceColor = [0.5 0.5 0.5];   % blue

        if(quantitative_data)
            b(2).FaceColor = [1 0.4 0];   % orange
        end

        % Label settings
        %set(gca, 'XTickLabel', electric_fields, 'XTickLabel', groups, 'XTickLabelRotation', 45);
        set(gca, 'XTickLabel', groups, 'XTickLabelRotation', 45, 'FontSize', 18, 'FontWeight', 'bold');

        ylabel('Directionality \Phi_x', 'FontSize',18);
        yticks([-1.0, -0.8, -0.6,-0.4, -0.2, 0])
        ylim([-1.0 0]);

        %legend({'Simulation', 'Experiment'}, 'Location', 'northwest');
        box on;
        set(gca,'fontsize',18);
        set(gcf, 'Color', 'white')
        exportgraphics(gcf, 'Directionality.pdf', 'ContentType', 'vector', 'Resolution',1000);
    end

    if(displacement3by2)
        tiledlayout(2,3);

        %plotting
        posxrun_str = ["posxr1", "posxr2", "posxr3", "posxr4", "posxr5", "posxr6"];
        posyrun_str = ["posyr1", "posyr2", "posyr3", "posyr4", "posyr5", "posyr6"];

        for i=1:6
            nexttile
            plot(cell_posData(1).(posxrun_str(i)), cell_posData(1).(posyrun_str(i)))
            hold on;
            xline(0, '-');
            yline(0, '-');
            xlabel('x-position (\mum)')
            ylabel('y-position (\mum)')
            %print1 = sprintf('EF Strength \n %.0f mV/mm', (fields(i)));
            %title(print1);
            xlim([-250, 150]);
            ylim([-150, 150]);
            set(gca,'fontsize',14);
            set(gcf, 'Color', 'white')
        end
        
        exportgraphics(gcf, 'sixPlots.pdf', 'ContentType', 'vector', 'Resolution',1000);    
    end
    
    num_density = 2;

    if(densityplotDIR)
        % compares directionality of different densities
        groups = {'15 mV/mm', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
        
        % Create grouped bar chart
        figure;
        
        b = bar(data);
        hold on

        % Set colors (optional)
        b(1).FaceColor = [0.5 0.5 0.5];     % blue
        b(2).FaceColor = [1 0.4 0];         % orange

        % Label settings
        set(gca, 'XTickLabel', groups, 'XTickLabelRotation', 45, 'FontSize', 18, 'FontWeight', 'bold');

        ylabel('Directionality \Phi_x', 'FontSize',18);
        yticks([-1.0, -0.8, -0.6,-0.4, -0.2, 0])
        ylim([-1.0 0]);

        legend({'Low Density', 'Higfh Density'}, 'Location', 'northwest');
        box on;
        set(gca,'fontsize',18);
        set(gcf, 'Color', 'white')

    end

    if(densityplotDISP)
        % compares displacement of different densities
        groups = {'15 mV/mm', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
        
        % Create grouped bar chart
        figure;
        
        b = bar(data);
        hold on

        % Set colors (optional)
        b(1).FaceColor = [0.5 0.5 0.5];     % blue
        b(2).FaceColor = [1 0.4 0];         % orange

        % Label settings
        set(gca, 'XTickLabel', groups, 'XTickLabelRotation', 45, 'FontSize', 18, 'FontWeight', 'bold');

        ylabel('Directionality \Phi_x', 'FontSize',18);
        yticks([-1.0, -0.8, -0.6,-0.4, -0.2, 0])
        ylim([-1.0 0]);

        legend({'Low Density', 'Higfh Density'}, 'Location', 'northwest');
        box on;
        set(gca,'fontsize',18);
        set(gcf, 'Color', 'white')
    end

end % end live conditional
end % end function