% Formal Plots for APS and Paper

%% NON-DENSITY; DC FIELD
% Required: [0, 30, 50, 75, 100, 200] V/m
% adjusting \mu tweaks calibration of disp/ dir
% CALIBRATING TO DIRECTIONALITY: \mu = 0.00065
close all;

% Directionality Histogram
groups = {'No EF', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
data_sim = zeros(length(fields),1);
err_sim = zeros(length(fields),1);

for i = 1:length(groups)
    data_sim(i,1) = cell_posData(i).directionality_mean;
    err_sim(i,1) = cell_posData(i).directionalitySD;
end

% Combine into matrix (columns = datasets, rows = groups)
data = data_sim;  % [6x1] matrix
err = err_sim;
% Create grouped bar chart
figure;
b = bar(data);

xBar = zeros(6,1);

hold on


xBar(:) = b(:).XEndPoints;


for i =1:length(groups)
    er = errorbar(xBar(i),data(i), err(i), 'o');
    er.Color = [0 0 0];
    hold on
end
hold on

% Set colors
b(1).FaceColor = [0.5 0.5 0.5];   % blue

% Label settings
set(gca, 'XTickLabel', groups, 'XTickLabelRotation', 45, 'FontSize', 18, 'FontWeight', 'bold');

ylabel('Directionality \Phi_x', 'FontSize',18);
yticks([-1.0, -0.8, -0.6,-0.4, -0.2, 0])
ylim([-1.0 0]);

leg_str = sprintf('Simulation, %d runs', runs);
%legend(leg_str, 'Location', 'southwest');
box on;
set(gca,'fontsize',18);
set(gcf, 'Color', 'white')


exportgraphics(gcf, 'Directionality_APS.pdf', 'ContentType', 'vector', 'Resolution',1000);





% Displacement Histogram
% CALIBRATING TO DISPLACEMENT:



groups = {'No EF', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
data_sim = zeros(length(fields),1);
err_sim = zeros(length(fields),1);

for i = 1:length(groups)
    data_sim(i,1) = cell_posData(i).displacement_mean;
    err_sim(i,1) = cell_posData(i).displacementSD;
end

% Combine into matrix (columns = datasets, rows = groups)
data = data_sim;  % [6x1] matrix
err = err_sim;
% Create grouped bar chart
figure;
b = bar(data);
hold on

xBar(:) = b(:).XEndPoints;

for i =1:length(groups)
    er = errorbar(xBar(i),data(i), err(i), 'o');
    er.Color = [0 0 0];
    hold on
end
hold on

% Set colors
b(1).FaceColor = [0.5 0.5 0.5];   % blue

% Label settings
set(gca, 'XTickLabel', groups, 'XTickLabelRotation', 45, 'FontSize', 18, 'FontWeight', 'bold');

ylabel('Displacement speed (\mum/min)', 'FontSize',15);
yticks([0.25, 0.5, 0.75, 1.0,1.25, 1.5,1.75, 2.0])
ylim([0.0 2.0]);

%leg_str = sprintf('Simulation, %d runs', runs);
%legend(leg_str, 'Location', 'northwest');
box on;
set(gca,'fontsize',15);
set(gcf, 'Color', 'white')

exportgraphics(gcf, 'Displacement_speed_APS.pdf', 'ContentType', 'vector', 'Resolution',1000);



%% NON-DENSITY AC
% REQUIRED: AC FIELD REVERSAL; NON-DENSITY

% choose 100 mV/mm case for analysis
cosTheta_200 = zeros(runTime, runs);
cosTheta_30 = zeros(runTime, runs);
cosTheta_100 = zeros(runTime, runs);
cosTheta_50 = zeros(runTime, runs);
cosTheta_75 = zeros(runTime, runs);

sumCellAnglex_200 = zeros(runTime, runs);
sumCellAnglex_100 = zeros(runTime, runs);
sumCellAnglex_30 = zeros(runTime, runs);
sumCellAnglex_50 = zeros(runTime, runs);
sumCellAnglex_75 = zeros(runTime, runs);

directionalityX_200 = zeros(runTime, runs);
directionalityX_100 = zeros(runTime, runs);
directionalityX_30 = zeros(runTime, runs);
directionalityX_50 = zeros(runTime, runs);
directionalityX_75 = zeros(runTime, runs);


for i = 1:runs
    cosTheta_200(:,i) = mean((cell_posData(i).direct6),2);
    cosTheta_100(:,i) = mean((cell_posData(i).direct5),2);
    cosTheta_30(:,i) = mean((cell_posData(i).direct2),2);
    cosTheta_50(:,i) = mean((cell_posData(i).direct3),2);
    cosTheta_75(:,i) = mean((cell_posData(i).direct4),2);


end

directionalityMean_200 = mean(cosTheta_200,2);
directionalityMean_100 = mean(cosTheta_100,2);
directionalityMean_30 = mean(cosTheta_30,2);
directionalityMean_50 = mean(cosTheta_50,2);
directionalityMean_75 = mean(cosTheta_75,2);


figure
for i = 1:runs
    s1 = scatter(2*time_control, cosTheta_30(:,i), "filled", 'd', "MarkerFaceColor", [0 0 1]);
    s1.MarkerFaceAlpha = 0.5;
    hold on

  %  s2 = scatter(2*time_control, cosTheta_50(:,i), "filled", 'd', "MarkerFaceColor", [0 0 1]);
   % s2.MarkerFaceAlpha = 0.5;
    hold on
    
 %   s3 = scatter(2*time_control, cosTheta_75(:,i), "filled", 'd', "MarkerFaceColor", [0 0 1]);
 %   s3.MarkerFaceAlpha = 0.5;
    hold on
    s4 = scatter(2*time_control, cosTheta_100(:,i), "filled", 'o', "MarkerFaceColor", [1 0 0]);
    s4.MarkerFaceAlpha = 0.5;
    hold on

    s5 = scatter(2*time_control, cosTheta_200(:,i), "filled", '^', "MarkerFaceColor", [1 0 1]);
    s5.MarkerFaceAlpha = 0.5;
    hold on
end

hold on
p5 = plot((time_control * 2), directionalityMean_200, '--', 'Linewidth', 3, 'Color', [0 0 0]);
hold on
p4 = plot((time_control * 2), directionalityMean_100, '-.', 'Linewidth', 3, 'Color', [0 0 0]);
hold on
p1 = plot((time_control * 2), directionalityMean_30, ':', 'Linewidth', 3, 'Color', [0 0 0]);
hold on
hold on
%p2 = plot((time_control * 2), directionalityMean_50, '-', 'Linewidth', 2, 'Color', [0.5 0 0.5]);hold on
hold on;
%p3 = plot((time_control * 2), directionalityMean_75, '--', 'Linewidth', 2, 'Color', [0 0.5 0.5]);
hold on;
xlabel('Time (minutes)', 'FontSize', 18);  ylabel('Directionality (\Phi_{x})', 'FontSize', 18);
%xticks = 2*(1:runTime)';
xticks([0 60 120 180 240])
xticklabels([0 60 120 180 240]);
xline((runTime / 2) * 2,'-.', 'FLIP', 'LineWidth',2, 'FontSize', 14);
ylim([-1.0,1.0]); xlim([0, runTime * 2]);

yline(0);
yline(1, '--');
yline(-1, '--');
txt1 = {'FIELD RIGHT'}; txt2 = {'FIELD LEFT'};
text((runTime*2 /9), 1.1, txt1, 'FontSize',14);
text((5*runTime /4), 1.1, txt2, 'FontSize',14);
set(gcf,'Color','white');
set(gca, 'FontSize', 18);

%legend([p1 p2 p3 p4 p5], {'30 V/m', '50 V/m', '75 V/m', '100 V/m', '200 V/m'}, 'Location','northwest')
legend([s1 s4 s5], {'30 V/m', '100 V/m', '200 V/m'}, 'Location','northwest')
exportgraphics(gcf, 'Field_flip_APS.pdf', 'ContentType', 'vector', 'Resolution',1000);


%% 6 plots for all EF strengths
% figure;
% tiledlayout(2,3);

%plotting
posxrun_str = ["posxr1", "posxr2", "posxr3", "posxr4", "posxr5", "posxr6"];
posyrun_str = ["posyr1", "posyr2", "posyr3", "posyr4", "posyr5", "posyr6"];

trajColorBlack = [ 0 0 0 ]; 
trajColorRed = [1 0 0];

Trajectory_data = struct();
trajectory_string = ["EF0", "EF30", "EF50", "EF75", "EF100", "EF200"];
Trajectory_red = struct();
Trajectory_black = struct();

for i = 1:6
    Trajectory_data.(trajectory_string(i)) = cell_posData(1).(posxrun_str(i));
    Trajectory_red.(trajectory_string(i)) = (Trajectory_data.(trajectory_string(i)) >0);
    Trajectory_black.(trajectory_string(i)) = (Trajectory_data.(trajectory_string(i)) <0);

    Trajectory_red.(trajectory_string(i)) = Trajectory_red.(trajectory_string(i))(end, :) .* ones(size(Trajectory_data.(trajectory_string(i))));
    Trajectory_black.(trajectory_string(i)) = Trajectory_black.(trajectory_string(i))(end, :) .* ones(size(Trajectory_data.(trajectory_string(i))));

end

for i=1:6
    nexttile

    % inner product with data to plot red
    Trajectory_data_tmp_x = Trajectory_red.(trajectory_string(i)) .* cell_posData(1).(posxrun_str(i));
    Trajectory_data_tmp_y = Trajectory_red.(trajectory_string(i)) .* cell_posData(1).(posyrun_str(i));

    plot(Trajectory_data_tmp_x, Trajectory_data_tmp_y, 'Color',trajColorRed, 'LineWidth', 0.7)
    hold on;

    % inner product with data to plot black
    Trajectory_data_tmp_x = Trajectory_black.(trajectory_string(i)) .* cell_posData(1).(posxrun_str(i));
    Trajectory_data_tmp_y = Trajectory_black.(trajectory_string(i)) .* cell_posData(1).(posyrun_str(i));

    plot(Trajectory_data_tmp_x, Trajectory_data_tmp_y, 'Color',trajColorBlack, 'LineWidth', 0.7)

    hold on;
    xline(0, '-');
    yline(0, '-');
    if(i == (4))
        xlabel('x-position (\mum)')
    end
    if(i == (5))
        xlabel('x-position (\mum)')
    end
    if(i == (6))
        xlabel('x-position (\mum)')
    end


    if(i == (1))
        ylabel('y-position (\mum)')
    end
    if(i == (4))
        ylabel('y-position (\mum)')
    end

print1 = title(sprintf('%.0f V/m', (fields(i))));

% Set the box and background properties
set(print1, 'Units', 'normalized', ...
            'Position', [0.8, 0.85, 0], ...
            'EdgeColor', 'black', ...      % Adds the box border
            'BackgroundColor', 'white', ... % Fills the box so plot lines don't cross text
            'LineWidth', 1, ...           % Thickness of the box border
            'Margin', 2);                 % Padding between text and the box border

    xlim([-250, 150]);
    ylim([-150, 150]);

    set(gca,'fontsize', 15);
    set(gcf, 'Color', 'white')
end




exportgraphics(gcf, 'sixPlots.pdf', 'ContentType', 'vector', 'Resolution',1000);


%% Density plots
if(density)
    % compares directionality of different densities
    groups = {'15 mV/mm', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
    dircrun_str = ["direct1", "direct2", "direct3", "direct4", "direct5", "direct6"];
    % Create grouped bar chart

    tmp_high = zeros(6, runs/2);
    tmp_low = zeros(6, runs/2);
    std_low = zeros(6, runs/2);
    std_high = zeros(6, runs/2);

    for i = 1:length(groups)
        for j = 1:runs
            if(mod(j,2))
                tmp_high(i, j) = cos(atan2((cell_posData(j).posymean(runTime, 1)), (cell_posData(j).posxmean(runTime, 1))));
            end
            if(~mod(j,2))
                tmp_low(i, j) = cos(atan2((cell_posData(j).posymean(runTime, 1)), (cell_posData(j).posxmean(runTime, 1))));
            end
        end
    end

    data = [low_density', high_density'];
    err = [LD_err', HD_err'];

    figure;

    b = bar(data);
    xBar = zeros(6,2);

    hold on
    for i = 1:2
        xBar(:, i) = b(i).XEndPoints;
    end

    for i =1:length(groups)
        for j = 1:2
            er = errorbar(xBar(i,j),data(i,j), err(i,j), 'o');
            er.Color = [0 0 0];
            hold on
        end
    end
    hold on

    % Set colors (optional)
    b(1).FaceColor = [0.5 0.5 0.5];     % blue
    b(2).FaceColor = [1 0.4 0];         % orange

    % Label settings
    set(gca, 'XTickLabel', groups, 'XTickLabelRotation', 45, 'FontSize', 18, 'FontWeight', 'bold');

    ylabel('Directionality \Phi_x', 'FontSize',18);
    yticks([-1.0, -0.8, -0.6,-0.4, -0.2, 0])
    ylim([-1.0 0]);

    legend({'Low Density', 'High Density'}, 'Location', 'southwest');
    box on;
    set(gca,'fontsize',18);
    set(gcf, 'Color', 'white')

end


%% single cell trajectory for supp.
single_cell = 0;
if(single_cell)

cell1_x  = cell_posData(1).posxr1(:,1); %#ok<UNRCH> 
cell1_y  = cell_posData(1).posyr1(:,1);
plot(cell1_x, cell1_y, 'color', [0 0 0], 'LineWidth', 3)
    hold on;
    xline(0, '-');
    yline(0, '-');
    xlabel('x-position (\mum)', 'FontSize',18)
    ylabel('y-position \mum)', 'FontSize',18)

    xlim([-250, 150]);
    ylim([-150, 150]);
end

%% gradient based trajectory plot 

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
exportgraphics(gcf, 'gradient.pdf', 'ContentType', 'vector', 'Resolution',1000);


