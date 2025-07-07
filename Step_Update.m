% Iteration update of parameters for function
function [xf, yf, vxf, vyf, cell_lifetime, tmp_time] = Step_Update(x0, y0, vx0, vy0, Fx, Fy, neibAngAvg, exempt, cell_lifetime, Cradius, growth_rate, tmp_time)

% Inputs: x(t),y(t) // vx(t),vy(t) // radius(t) // Fx(t),Fy(t) // alignment angle(t)
% Outputs: x(t+1), y(t+1) // vx(t+1), vy(t+1) // radius(t+1)

% Declaration of Globals
global NumCells dt eta vels_med critRad runTime speed_decay dim2noise dim1noise etaX etaY                                      %#ok<GVMIS> 

%% Preallocation of updated values
vxf = zeros(NumCells,1);
vyf = zeros(NumCells,1);
vxNat = zeros(NumCells, 1);
vyNat = zeros(NumCells, 1);
angNatural = zeros(NumCells, 1);
xf = x0;
yf = y0;

% Loop for updating the values of cell
for i = 1:NumCells
    %%  Position Update

    % Considers all cells still alive 
    if(exempt(i,1))
        cell_lifetime(i,1) = cell_lifetime(i,1) + 1;
        xf(i,1) = x0(i,1) + vx0(i,1) *dt;
        yf(i,1) = y0(i,1) + vy0(i,1)*dt;

        %% Velocity update
        % Natural angle: mean angle of self and friends, influenced by noise
        if(dim1noise)
        angNatural(i,1) = neibAngAvg(i,1) + (eta * (rand() - 0.5)*pi)*sqrt(dt);

        % natural velocity of cell (excludes CC & CF)
        vxNat(i,1) = vels_med * cos(angNatural(i,1));
        vyNat(i,1) = vels_med * sin(angNatural(i,1));
        end

        if(dim2noise)
            angNatural(i,1) = neibAngAvg(i,1);
            vxNat(i,1) = vels_med * cos(angNatural(i,1) + etaX * (rand() - 0.5)*pi*sqrt(dt));
            vyNat(i,1) = vels_med * sin(angNatural(i,1) + etaY * (rand() - 0.5)*pi*sqrt(dt));
        end

        % New velocity vector based on how interaction forces affected angles
        % (componentwise)
        vxf(i, 1) = (vxNat(i,1) + Fx(i,1))* dt; % x component
        vyf(i, 1) = (vyNat(i,1) + Fy(i,1))* dt;% y component
        if(Cradius(i,1) > 0.9 * critRad)
            tmp_time(i,1) = tmp_time(i,1) + 1;
            pow = -(tmp_time(i,1) / (speed_decay)); 
            vxf(i,1) = vxf(i,1) * exp(pow);
            vyf(i,1) = vyf(i,1) * exp(pow);
        end
    end % end live cell condition
end % end all cells loop
end % end function



