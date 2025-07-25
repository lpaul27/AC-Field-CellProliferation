function            [Cradius,x, y, vx, vy, vel_ang, x_time, y_time, theta_time, RadTracker, R, G, B, Pressure, exempt, vx_time, vy_time, cell_lifetime, tmp_time]...
    = RadGrowth(Cradius0, Pressure, x, y, vel_ang, vx, vy, x_time, y_time, time, theta_time, RadTracker, R, G, B, exempt, vx_time, vy_time, cell_lifetime, Ex_strength, Ey_strength, growth_rate, tmp_time)
%% Allows cell growth, mitosis and death

% declaration of constants
global dt NumCells critRad critical_pressure vels_med daughter_noise ...
    death_pressure death_rate chill runTime Field rand_division                                                              %#ok<GVMIS>

% Initialization
Cradius = Cradius0;
death_chance = zeros(NumCells, 1);
field_angle = atan2(Ey_strength, Ex_strength);

%% loop over all cells
for i = 1:NumCells
    %% Growth conditional

    if(Pressure(i,1) < critical_pressure && Cradius(i,1) <= critRad && exempt(i,1))
        Cradius(i,1) = Cradius0(i,1) + (growth_rate(i, 1) ./ (2*pi.*Cradius0(i,1))) * dt;
    end % end growth conditional

    %% Mitosis Conditional

    if(Cradius(i,1) >= critRad && Pressure(i,1) <= critical_pressure && exempt(i,1))
        % Add parameters if mitosis criterion is reached
        NumCells = NumCells + 1;
        % random based position mitosis
        if(~Field || rand_division)
            x(NumCells,1) = x(i,1) + Cradius(i,1)*(1 - 1/sqrt(2)) * (rand() - 0.5);
            y(NumCells,1) = y(i,1) + Cradius(i,1)*(1 - 1/sqrt(2)) * (rand() - 0.5);
            x(i,1) = x(i,1) + Cradius(i,1)*(1 - 1/sqrt(2)) * (rand() - 0.5);
            y(i,1) = y(i,1) + Cradius(i,1)*(1 - 1/sqrt(2)) * (rand() - 0.5);    
            y(NumCells,1) = y(i,1) + Cradius(i,1)*(1 - 1/sqrt(2)) * (rand() - 0.5);
            vx(NumCells, 1) = vels_med * cos(vel_ang(i) + daughter_noise * pi * (rand() - 0.5));
            vy(NumCells, 1) = vels_med * sin(vel_ang(i) + daughter_noise * pi * (rand() - 0.5));
        end
        % perpendicular to field based mitosis
        if(Field && ~rand_division)
            x(NumCells, 1) = x(i,1) + (1 - 1/sqrt(2))*Cradius(i,1) * cos(field_angle * daughter_noise * (rand() - 0.5));
            y(NumCells, 1) = y(i,1) - (1 - 1/sqrt(2))*Cradius(i,1) * sin(field_angle * daughter_noise * (rand() - 0.5));
            x(i,1) = x(i,1) - (1 - 1/sqrt(2))*Cradius(i,1)*rand() * cos(field_angle * daughter_noise * (rand() - 0.5));
            y(i, 1) = y(i,1) + (1 - 1/sqrt(2))*Cradius(i,1)*rand() * sin(field_angle * daughter_noise * (rand() - 0.5));
            
            vx(NumCells, 1) = vels_med * cos(field_angle) + daughter_noise * pi * (rand() - 0.5);
            vy(NumCells, 1) = vels_med * sin(field_angle) + daughter_noise * pi * (rand() - 0.5);
            vx(i,1) = -vels_med * cos(field_angle) + daughter_noise * pi * (rand() - 0.5);
            vy(i,1) = -vels_med * sin(field_angle) + daughter_noise * pi * (rand() - 0.5);
        end
        
        vel_ang(NumCells, 1) = atan2(vy(NumCells,1), vx(NumCells,1));
        Cradius(NumCells, 1) = Cradius(i,1) ./ sqrt(2);
        Cradius(i,1) = Cradius(i,1) ./ sqrt(2);
        Pressure(NumCells, 1) = 0;
        exempt(NumCells , 1) = 1;

        % Update Plotting parameters
        R(NumCells,1) = 0;
        G(NumCells,1) = 0.7;
        B(NumCells, 1) = 1;

        % Update tracker to account for new cells
        x_time(1:time, NumCells) = x_time(1:time, i);  
        y_time(1:time, NumCells) = y_time(1:time, i);
        theta_time(1:time, NumCells) = theta_time(1:time, i);
        RadTracker(1:time, NumCells) = RadTracker(1:time, i);
        vx_time(1:time, NumCells) = vx_time(1:time, i);
        vy_time(1:time, NumCells) = vy_time(1:time, i);
        cell_lifetime(NumCells, 1) = 0;
        cell_lifetime(i, 1) = 0;
        tmp_time(i,1) = 0;
        tmp_time(NumCells,1) = 0;

    end % end mitosis conditional

    %% Death conditional

    % apoptosis random death based on length of cell life
    death_chance(i,1) = (death_rate / runTime) .*cell_lifetime(i,1);
    tmp = rand();
    if((Pressure(i,1) >= death_pressure || tmp < death_chance(i, 1)) && exempt(i,1) && time > chill)

        % Mark cell as exempt from simulation (dead)
        exempt(i,1) = 0;

        % update parameters
        vx(i, 1) = 0;
        vy(i, 1) = 0;
        vel_ang(i, 1) = atan2(vy(i,1), vx(i,1));
        Pressure(i,1) = 0;

        % Set color to pure red for a dead cell
        R(i,1) = 0;
        G(i,1) = 0;
        B(i, 1) = 0;
    end % end cell death region
    %% Update RGB scale
    % cells which are still alive are colored based on amount of pressure
    if(exempt(i,1))
        B(i,1) = 1 - (Pressure(i,1) ./ (critical_pressure + Pressure(i,1)));
        G(i,1) = B(i,1);
        R(i,1) = 1 - B(i,1);
    end % end color conditional
end % end cell loop
end % end function