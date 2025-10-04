%% Main file for cell movement in an AC field
% Conatins declaration, parameters, time loop, all function calls

clear all;                                                                  %#ok<CLALL>
close all;
%% Parameters for model
% Global parameters declaration
global NumCells dt lbox vels_med eta nu neighborWeight k R_boundary Cell_radius ...
    c_rec c_lig adh runTime Field xphi yphi w ExMax EyMax mu ...
    critRad Ccyclet critical_pressure daughter_noise Cell_std death_rate ...
    death_pressure chill dim2directionality displacement live polarhist dim1directionality...
    Discrete Sine dim1displacement rand_division speed_decay dim1noise dim2noise etaX etaY ...                                                                                        %#ok<GVMIS>
    disphist velocity_noise sigmax sigmay directednessplot velocity_mag_noise displacement3by2 ...                                                                                               %#ok<GVMIS>
    densityplotDIR densityplotDISP densityDirTime tau E0

%% Types of plots
live = 0;                               % Enables Live Visualization
dim1directionality = 0;                 % enables 1D Directionality plot
dim2directionality = 0;                 % Enables 2D Directionality plot
dim1displacement = 0;                   % Enables 1D Displacement plot
displacement = 0;                       % Enables 2D Displacement plot
polarhist = 0;                          % Enables polar histogram plot
disphist =1;                           % enables displacement histogram plot
directednessplot = 0;                   % enables directedness histogram plot
displacement3by2 = 1;                   % enables displacement of all fields in a 3x2
densityplotDIR = 0;                     % enables density plot for directionality
densityplotDISP = 0;                   % enables density plot for displacement
densityDirTime = 0;                     % enables directionality over time plot for density

% Number of runs to be averaged across
runs = 3;

% begin start timer
tStart = tic;

%% constants found to make displacement plost match
noises = [0.26, 0.24, 0.2, 0.28, 0.32, 0.56];
forceScale = [1,1.11,1.4,0.85,0.77,0.45];
noiseX = [0.28,0.28,0.23,0.33,0.37,0.61];
v0 = [0.083, 0.11, 0.185, 0.07, 0.059, 0.026];
nufric = [1, 1.15, 1.43, 0.93, 0.81,0.46];

%% Declaration of strings for class for access in loops
posxrun_str = ["posxr1", "posxr2", "posxr3", "posxr4", "posxr5", "posxr6"];
posyrun_str = ["posyr1", "posyr2", "posyr3", "posyr4", "posyr5", "posyr6"];
Directrun_str = ["direct1", "direct2", "direct3", "direct4", "direct5", "direct6"];

%% Outer loop to control field strengths:
% fields (normal)           densityfields(density):
% 1: NO EF                  15 V/m
% 2: 30 V/m                 30 V/m
% 3: 50 V/m                 50 V/m
% 4: 75 V/m                 75 V/m
% 5: 100 V/m                100 V/m
% 6: 200 V/m                200 V/m

fields = [0, 30, 50, 75,  100, 200];
densityfields = [15, 30, 50, 75, 100, 200];
density = false;        % signals that a varying number of cells is intended
% *default is 50 cells 

%% Preallocation
directedness = zeros(length(fields),1);
dispAvg = zeros(length(fields),1);
cell_posData = struct();
cell_densitydata = struct();

% Set to a higher value if multiple field intended
for z = 1:6
    % Preallocation for data storage
    displacement_run = zeros(runs, 1);
    displacementAvg = zeros(runs,  1);
    dispruntheta = zeros(runs,1);

    %% Inner loop for number of runs at desired field strength
    for p = 1:runs
        % Begin Simulation timer
        %% Domain Parameters
        runTime = 120;                           % total runTime of simulation
        dt = 1;                                  % time step
        NumCells = 50;
        if(density)
            if(~mod(p,2))
                NumCells =35;                           % number of cells in simulation
            end
            if(mod(p,2))
                NumCells = 210;
            end
        end
        vels_med = 0.7;                        % initial velocity param center point
        lbox = 1050;                             % size of the box particles are confined to
        R_boundary = 150; %lbox/6;                     % Sample domain size for cells to begin
        chill = 15;                              % chill time to suppress cell death

        %% Cell Parameters
        critRad = 12;                            % critical radius for mitosis
        Ccyclet = 1000;                      % benchmark cell cycle time
        death_rate = 1e-200;                     % Cell death rate
        death_pressure = 1000;                   % Pressure required for apoptosis
        critical_pressure = 0.05;                % Critical presssure for dormancy
        Cell_radius = 10;                        % fixed cell radius
        Cell_std = 0.08;                         % Standard Deviation of cell radii
        speed_decay = 100;                       % speed decay rate for mitosis

        %% Cell-cell parameters
        k = 0.01;                               % constant in force repulsion calculation (~elasticity)
        daughter_noise = 0.1;                   % noise strength in mitosis separation
        nu = 1;                                 % friction factor
        mu = 1;                                 % electrical mobility
        neighborWeight = 1;                  % group movement weighting
        c_rec = 0.9;                            % mean receptor concentration (normalized)
        c_lig = 0.9;                            % mean ligand concentration (normalized)
        adh = 0;                                % adhesive coefficient

        %% Cell-Field parameters
        % Discrete Parameters
        Field = 1;                              % Signals to time varying fields that field is on if 1
        rand_division = 0;                      % Enables field-directed mitosis
        Discrete = 1;                           % Enables Discrete field change

        ExMax = fields(z) / 2500;               % x field max
        EyMax = 0;                              % y field max
        absE = sqrt(EyMax^2 + ExMax^2);         % magnitude of field
        vels_med = vels_med*(1+fields(z)/200);  % initial velocity param center point


        % Sinusoidal parameters
        % f(t) = A sin(wt + o)                  % form
        Sine = 0;                               % enables sinusoidal field change
        w = 8*pi /(runTime);                    % angular frequency
        xphi = 0;                               % x field offset
        yphi = 0;                               % y field offset

        %% Simulation type parameters
        dim1noise = 1;                          % signals type of noise (1D)
        %eta = noise * (1+ alpha / (emax/absE + 1));
        eta = 0.15;
        dim2noise = 0;                          % signals type of noise (2D)
            etaX = eta / 2;                     % X component of noise strength
            etaY = 3*eta;                       % Y component of noise strength
            velocity_noise = 0;                 % signals type of noise (velocity)
            sigmax = eta * ExMax / (1 + absE);  % velocity based noise parameter (x)
            sigmay = eta * EyMax / (1 + absE);  % velocity based noise parameter (y)
            velocity_mag_noise = 0;             % signals type of noise (velocity magnitude)
        
        % Polarity (intracellular)
        intracellular_signalling = 1;
            tau = 20;                            % time scale parameter
            E0 = 50;                           % Relative reference point for saturation
            

        %% Plot parameters

        time_control = (1:runTime)';            % time axis for plotting
        R = zeros(NumCells, 1);                 % Red scale for plotting
        G = ones(NumCells, 1);                  % Green scale for plotting
        B = ones(NumCells, 1);                  % Blue scale for plotting



        %% Initialization of Variables
        % Preallocates values for optimal computation

        x_time = zeros(runTime, NumCells);      % Matrix of x position for each step
        y_time = zeros(runTime, NumCells);      % Matrix of y position for each step
        vx_time = zeros(runTime, NumCells);     % Matrix of vx position for each step
        vy_time = zeros(runTime, NumCells);     % Matrix of vy position for each step
        theta_time = zeros(runTime, NumCells);  % Matrix of angle for each step
        timer = zeros(runTime, 3);              % Timer to keep track of computational efficiency
        RadTracker = zeros(runTime, NumCells);  % tracker of cell size
        exempt = ones(NumCells, 1);             % cell death logical
        cell_lifetime = zeros(NumCells, 1);     % cell lifetime tracker
        growth_rate = zeros(NumCells, 1);       % initializes growth rate for all cells
        tmp_time = zeros(NumCells, 1);
        EF_x = zeros(NumCells, 1);
        EF_y = zeros(NumCells, 1);
        Epressure = zeros(NumCells, 1);

        %% Plotting Parameters
        % Parameters for live simulation visualization

        if(live)
            cell=figure;                                                            %#ok<UNRCH>
            cell.WindowState = 'maximized';
            axis([0 lbox 0 lbox])
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',12);
            axis('square')
            hold on
            skip_points = 14;
        end
        %% Initialization of System
        % Based on Monte Carlo initialization

        % Initialize cells
        [x, y, vx, vy, Cradius, vel_ang, polarity] = Initialize();

        %% Simulation loop
        for time = 1:runTime
            % Time loop
            if(exempt == zeros(size(exempt)))
                % if all cells die, end simulation
                disp('All cells dead')
                break
            end
            for j = 1:NumCells
                growth_rate(j, 1) = (pi * randgaussrad(critRad, (critRad / 2)).^2) / (2* Ccyclet);
            end
            % Reinitializes electric field
            if(time || Field)
                [u, v, X, Y, Ex_strength, Ey_strength] = EF_Grid_Init(time);
            end
            % Store plotting data

            % Displacement data
            x_time(time, :) = x(:,1);
            y_time(time, :) = y(:,1);

            % Velocity data
            vx_time(time, :) = vx(:, 1);
            vy_time(time,:) = vy(:,1);

            % Directionality and Mitosis data
            theta_time(time, :) = vel_ang(:,1);
            RadTracker(time, :) = Cradius(:,1);


            %% Force update functions (cell-cell & cell-field)

            % begin cell-cell timer
            CCtimer = tic;

            % cell-cell force function
            [Fx, Fy, neibAngAvg, Cpressure] = Interaction_Forces...
                (x, y, Cradius, vel_ang, exempt);

            timer(time, 1) = toc(CCtimer);
            % end cell-cell timer



            % begin cell-field timer
            CFtimer = tic;

            % cell-field force function
            if(~intracellular_signalling)
                [EF_x, EF_y, Epressure] = Electric_Force(Cradius, x, y, u, v, X, Y);
            end
            if(intracellular_signalling)
                [polarity, vx, vy] = Polarization(polarity, ExMax, EyMax, vel_ang, time);
            end

            timer(time, 2) = toc(CFtimer);
            % end cell-field timer

            %% Step update

            % Begin step update timer
            Steptimer = tic;

            % Calculate net force; exclude dead cells
            Fx_net = nu*(Fx + mu*EF_x) .* exempt;
            Fy_net = nu*(Fy + mu*EF_y) .* exempt;

            % Calculate the net pressure
            Pressure = Epressure + Cpressure;

            % Step update function
            [x, y, vx, vy, cell_lifetime, tmp_time] = Step_Update(x, y, vx, vy, Fx_net, Fy_net, neibAngAvg, exempt, cell_lifetime, Cradius, growth_rate, tmp_time);
            vel_ang = atan2(vy,vx);

            % end step update timer
            timer(time,3) = toc(Steptimer);

            % Update cell size // Mitosis // Apoptosis
            [Cradius,x, y, vx, vy, vel_ang, x_time, y_time, theta_time, RadTracker, R, G, B, Pressure, exempt, vx_time, vy_time, cell_lifetime, tmp_time, polarity]...
                = RadGrowth...
                (Cradius, Pressure, x, y, vel_ang, vx, vy, x_time, y_time, time, theta_time,RadTracker, R, G, B, exempt, vx_time, vy_time, cell_lifetime, Ex_strength, Ey_strength, growth_rate, tmp_time, polarity);

            %% Live Simulation visualization plot
            if(live)
                scale_efield = 2;                                                   %#ok<UNRCH>
                x_efield_plot = reshape(X,length(X)^2,1);
                y_efield_plot = reshape(Y,length(Y)^2,1);
                u_efield_plot = reshape(u,length(u)^2,1);
                v_efield_plot = reshape(v,length(v)^2,1);
                v_result = [vx vy];
                v_result_norm = sqrt(diag(v_result * v_result'));

                cla
                set(gcf,'doublebuffer','on')
                hold on;
                skip_nth =51;
                quiver(x_efield_plot(1:skip_nth:end),y_efield_plot(1:skip_nth:end),scale_efield*u_efield_plot(1:skip_nth:end),scale_efield*v_efield_plot(1:skip_nth:end), 'Color', [1, 0., 0],   'LineWidth', 1., 'MaxHeadSize', 0.9);
                hold on;
                quiver(x,y,Cradius .* vx./(0.5*v_result_norm),Cradius.*vy./(0.5*v_result_norm), 'Color',[0, 0, 0], 'MarkerSize', 10, 'LineWidth', 1.5,  'AutoScale', 'off') ;
                hold on;
                for i = 1:NumCells
                    circles(x(i), y(i), Cradius(i), 'facecolor', [R(i), G(i), B(i)]);
                end
                hold on;
                drawnow
                hold on
            end
        end % end time loop
        %% Static Plot Visualization
        x_centered = abs(x_time - x_time(1,:));
        y_centered = abs(y_time - y_time(1,:));
        displacement_mag = sqrt(x_centered.^2+y_centered.^2);
        %         displacementrunx = (xavg(runTime, 1) - xavg(1,1)) / runTime;
        %         displacementruny = (yavg(runTime, 1) - yavg(1,1)) / runTime;
        displacement_run(p,1) = mean(displacement_mag(runTime, :),2) / runTime;

        x_raw = (x_time - x_time(1,:));
        y_raw = (y_time - y_time(1,:));
        dispruntheta = mean(cos(atan2(y_raw(runTime,:) , x_raw(runTime, :))),2);
        posxrun_str = ["posxr1", "posxr2", "posxr3", "posxr4", "posxr5", "posxr6"];
        posyrun_str = ["posyr1", "posyr2", "posyr3", "posyr4", "posyr5", "posyr6"];

        %z = 1;              % resets field strength to index correctly
        %% Important: Leave z as the index of the loop if doing all field strengths, otherwise, 
        % it must be indexed to 1

        switch z
            case 1
                cell_posData(p).posxr1 = x_time - x_time(1,:);
                cell_posData(p).posyr1 = y_time - y_time(1,:);
                cell_posData(p).theta1 = theta_time;
                y_tmp = cell_posData(p).posyr1;
                x_tmp = cell_posData(p).posxr1;
                cell_posData(p).direct1 = cos(atan2(y_tmp, x_tmp));
            case 2
                cell_posData(p).posxr2 = x_time - x_time(1,:);
                cell_posData(p).posyr2 = y_time - y_time(1,:);
                cell_posData(p).theta2 = theta_time;
                y_tmp = cell_posData(p).posyr2;
                x_tmp = cell_posData(p).posxr2;
                cell_posData(p).direct2 = cos(atan2(y_tmp, x_tmp));
            case 3
                cell_posData(p).posxr3 = x_time - x_time(1,:);
                cell_posData(p).posyr3 = y_time - y_time(1,:);
                cell_posData(p).theta3 = theta_time;
                y_tmp = cell_posData(p).posyr3;
                x_tmp = cell_posData(p).posxr3;
                cell_posData(p).direct3 = cos(atan2(y_tmp, x_tmp));
            case 4
                cell_posData(p).posxr4 = x_time - x_time(1,:);
                cell_posData(p).posyr4 = y_time - y_time(1,:);
                cell_posData(p).theta4 = theta_time;
                y_tmp = cell_posData(p).posyr4;
                x_tmp = cell_posData(p).posxr4;
                cell_posData(p).direct4 = cos(atan2(y_tmp, x_tmp));
            case 5
                cell_posData(p).posxr5 = x_time - x_time(1,:);
                cell_posData(p).posyr5 = y_time - y_time(1,:);
                cell_posData(p).theta5 = theta_time;
                y_tmp = cell_posData(p).posyr5;
                x_tmp = cell_posData(p).posxr5;
                cell_posData(p).direct5 = cos(atan2(y_tmp, x_tmp));
            case 6
                cell_posData(p).posxr6 = x_time - x_time(1,:);
                cell_posData(p).posyr6 = y_time - y_time(1,:);
                cell_posData(p).theta6 = theta_time;
                y_tmp = cell_posData(p).posyr6;
                x_tmp = cell_posData(p).posxr6;
                cell_posData(p).direct6 = cos(atan2(y_tmp, x_tmp));
        end


    end
    dispAvg(z,1) = mean(displacement_run);
    directedness(z,1) = (dispruntheta);
end


for i = 1:length(1)
    for j = 1:runs
        cell_posData(i).posxmean = mean(cell_posData(j).(posxrun_str(i)),2);
        cell_posData(i).posxSD = std(cell_posData(j).(posxrun_str(i)),0,2);

        cell_posData(i).posymean = mean(cell_posData(j).(posyrun_str(i)),2);
        cell_posData(i).posySD = std(cell_posData(j).(posyrun_str(i)),0,2);

        tmp = cell_posData(j).(Directrun_str(i));
        cell_posData(i).directionality_mean = mean(tmp(runTime,:),2);
        cell_posData(i).directionalitySD = std(tmp(runTime,:),0,2) / sqrt(runs);
    end
end

Visualize(x_time,y_time, theta_time, time_control, dispAvg, directedness, fields, cell_posData);

% end timer of full sequence
toc(tStart);
% end simulation