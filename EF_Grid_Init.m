% Void Function to initialize Electric Field grid
% Allows a discretized creation of the field at equidistant intervals in a 2x2 grid

% Inputs: Domain size, time, total time of simulation
% Outputs: Electric field vector at each 'integer' coordinate (m,n), grid
% coordinates of field

function [u, v, X, Y, Ex_strength, Ey_strength] = EF_Grid_Init(time)

global lbox runTime Field xphi yphi w ExMax EyMax Discrete Sine                                  %#ok<*NUSED,GVMIS>
% Electric field represented as a function of position

%% Case 1: Representation based on experimentalist's data:
% Uniform Electric Field and discontinuous change to opposite polarity
% at uniform time intervals
if(Discrete && Field)
    if(time > (runTime / 2))
        Ex_strength = -ExMax;
        Ey_strength = -EyMax;
    else
        Ex_strength = ExMax;
        Ey_strength = EyMax;
    end % end field change cond
else
    Ex_strength = ExMax;
    Ey_strength = EyMax;
end % end discrete conditional
%% Case 2: Representation based on a sine function
if(Sine && Field)
    % Time varying magnitude
    Ex_strength = ExMax * sin((w * time) + xphi);
    Ey_strength = EyMax * sin((w*time) + yphi);
end % end sine cond

%% Update grid
% Field off
if(~Field && ~Discrete)
    Ex_strength = 0;
    Ey_strength = 0;
end

[X,Y] = meshgrid(0.1:2:lbox+0.5);
u=(Ex_strength)*ones(size(X));
v= (Ey_strength)* ones(size(Y));
end % end function