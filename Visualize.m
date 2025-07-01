function [] = Visualize(x_time,y_time, theta_time, time_control)
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
global NumCells runTime displacement live dim2directionality polarhist dim1directionality dim1displacement         %#ok<GVMIS>

% runs nothing if live simulation is on
if(~live)
    %% Displacement track graph
    % Visualization of cell proliferation

    if(displacement)

        %calculations
        x_avg = mean(x_time, 2);
        y_avg = mean(y_time, 2);

        %plotting
        figure
        hold on;
        plot((x_time - x_time(1,:)), (y_time - y_time(1,:)))
        hold on;
        plot((x_avg - x_avg(1,1)), (y_avg - y_avg(1,1)), 'w', 'LineWidth', 3);
        xline(0, '-');
        yline(0, '-');
        xlabel('x-displacement (a.u)')
        ylabel('y-displacement (a.u)')
    end

    %% 1D Directionality Graph
    % Visualization of allignment to 180 degree field change
    if(dim1directionality)
        cosTheta = cos(theta_time);
        sumCellAnglex = sum(cosTheta,2);
        directionalityX = sumCellAnglex ./ NumCells;
        confIntXrise = find(directionalityX >= (max(directionalityX) - 0.05));
        confIntXfall = find(directionalityX <= (min(directionalityX(10:end)) + 0.05) & gradient(directionalityX) < 0);
        XriseTime = confIntXrise(1,1);
        XfallTime = confIntXfall(1,1);

        figure
        plot(time_control, directionalityX)
        xlabel('Time (steps)');  ylabel('Directionality (\Phi)');
        xline((runTime / 2),'-.', 'TURN', 'LineWidth',2);
        ylim([-1.2,1.2]); xlim([0, runTime]);
        print1 = sprintf('X Rise: %f', XriseTime);
        print2 = sprintf('X Fall: %f', (XfallTime - (runTime / 2)));
        xline(XriseTime, ':', print1, 'Color', 'b', 'LineWidth',1);
        xline(XfallTime, ':', print2, 'Color', 'b', 'LineWidth',1);
        yline(0);
        txt1 = {'FIELD RIGHT'}; txt2 = {'FIELD LEFT'};
        text((runTime /4), 1.1, txt1);
        text((5*runTime /8), 1.1, txt2);

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
        title('Field Direction: +y')
    end % end single axis directionality
end % end live conditional
end % end function