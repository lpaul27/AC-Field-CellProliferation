% compares directionality of different densities
groups = {'15 mV/mm', '30 mV/mm', '50 mV/mm', '75 mV/mm', '100 mV/mm', '200 mV/mm'};
dircrun_str = ["direct1", "direct2", "direct3", "direct4", "direct5", "direct6"];
posxrun_str = ["posxr1", "posxr2", "posxr3", "posxr4", "posxr5", "posxr6"];
posyrun_str = ["posyr1", "posyr2", "posyr3", "posyr4", "posyr5", "posyr6"];
% Create grouped bar chart

cell_vector_lowDensity = struct('xcomponent', 'ycomponent');
cell_vector_highDensity = struct('xcomponent', 'ycomponent');

Phi_rho = zeros(6,10);
Phi_rho_low = zeros(6,5);
Phi_rho_high = zeros(6,5);
Phi_rho_SD = zeros(6,10);


for i = 1:6
    for j = 1:10
        tmpx = cell_posData(j).(posxrun_str(i))(120, :);
        tmpy = cell_posData(j).(posyrun_str(i))(120, :);

        tmpatan2 = atan2(tmpy, tmpx);
        tmp_Phi = cos(tmpatan2);

        tmp_avgPhi = mean(tmp_Phi);

        Phi_rho(i,j) = tmp_avgPhi;
        Phi_rho_SD(i,j) = std(tmp_Phi);
        
        if(mod(j, 2))
            Phi_rho_high(i,round(j/2)) = Phi_rho(i,j);
        end
        if(~mod(j,2))
            Phi_rho_low(i,round(j/2)) = Phi_rho(i,j);
        end
    end
end
disp(tmpx)
Phi_lowDensity = mean(Phi_rho_low, 2);
Phi_lowDensity_SD = std(Phi_rho_low,0, 2);

Phi_highDensity = mean(Phi_rho_high, 2);
Phi_highDensity_SD = std(Phi_rho_high,0, 2);

data = [Phi_highDensity, Phi_lowDensity];

figure;

b = bar(data);
xBar = zeros(6,2);
err = zeros(6,2);
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

legend({'High Density', 'Low Density'}, 'Location', 'southwest');
box on;
set(gca,'fontsize',18);
set(gcf, 'Color', 'white')
