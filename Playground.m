figure;
tiledlayout(1,2);

%plotting
posxrun_str = ["posxr1", "posxr6"];
posyrun_str = ["posyr1", "posyr6"];
fields = [0, 200];
trajectory_string = ["EF0", "EF30", "EF50", "EF75", "EF100", "EF200"];

t = (linspace(0, 1, runTime))';
cmap = zeros(runTime-1, 3);
for i=1:2
    nexttile

    % inner product with data to plot red
    Trajectory_data_tmp_x=cell_posData(1).(posxrun_str(i));
    Trajectory_data_tmp_y = cell_posData(1).(posyrun_str(i));


    for j = 1:(runTime-1)
        for k = 1:10
            R = 0.3 + 0.5*t(j);
            G = 0.7 * t(j);
            B = 0.5 - 0.3* t(j);
            cmap(j, :) = [R G B];

            plot([Trajectory_data_tmp_x(j, k) Trajectory_data_tmp_x(j+1, k)], [Trajectory_data_tmp_y(j, k)  Trajectory_data_tmp_y(j+1, k)], 'Color', [R G B], 'LineWidth', 2.0)
            hold on;
        end
    end 

    hold on;
    xline(0, '-');
    yline(0, '-');
    if(i == (1))
        xlabel('x-position (\mum)')
    end
    if(i == (2))
        xlabel('x-position (\mum)')

        colormap(cmap);
        clim([0 3])
        cb = colorbar;
        cb.Label.String = 'Time (hours)';
        cb.Ticks = [0 1 2 3];
    end
    if(i == (1))
        ylabel('y-position (\mum)')
    end



print1 = title(sprintf('%.0f V/m', (fields(i))));

% Set the box and background properties
set(print1, 'Units', 'normalized', ...
            'Position', [0.9, 0.95, 0], ...
            'EdgeColor', 'black', ...      % Adds the box border
            'BackgroundColor', 'white', ... % Fills the box so plot lines don't cross text
            'LineWidth', 1, ...           % Thickness of the box border
            'Margin', 2);                 % Padding between text and the box border

    xlim([-250, 150]);
    ylim([-150, 150]);

    set(gca,'fontsize', 15);
    set(gcf, 'Color', 'white')
end
fields = [0, 30, 50, 75,  100, 200];

