% Initialization of system void
% Calls for initialization of system

% Inputs: Domain size, avg velocity, total number of cells
% Outputs: position of each cell, component velocity of each cell, radius
% of each cell, angle of travel

function [x, y, vx, vy, Cradius, vel_ang, polarity] = Initialize()

% Declaration of constants
global lbox vels_med NumCells R_boundary Cell_radius Cell_std               %#ok<GVMIS>

%% Initialization of variables
Cradius = zeros(NumCells, 1);
x = zeros(NumCells, 1);
y = zeros(NumCells, 1);
vx = zeros(NumCells, 1);
vy = zeros(NumCells, 1);
r = zeros(NumCells, 1);
theta = zeros(NumCells, 1);
vel_ang = zeros(NumCells, 1);
speed = zeros(NumCells, 1);
polarity = zeros(NumCells, 1);

% Preallocate space
maxAttempts = 10000;      % Safety cap for how many random points to test
centers = NaN(NumCells, 2);  % Stores accepted circle centers
theta = 2*pi*rand(maxAttempts, 1);

% Generate a large batch of random points
% X= R_boundary*sqrt(rand(maxAttempts,1)).* cos(theta)+lbox/2;
% Y= R_boundary*sqrt(rand(maxAttempts,1)).* sin(theta)+lbox/2;
% points = [X Y];

% First point is always accepted
% accepted = points(1, :);
% acceptedCount = 1;

% Vectorized filtering of remaining points
% for i = 2:maxAttempts
%     dists = sqrt(sum((accepted - points(i,:)).^2, 2));
%     if all(dists >= 2*Cell_radius)
%         accepted = [accepted; points(i,:)];
%         acceptedCount = acceptedCount + 1;
%         if acceptedCount == NumCells
%             break
%         end
%     end
% end
% x = accepted(:,1);
% y = accepted(:,2);

%% Iteration loop over all cells to initialize position and velocity
for cells = 1:NumCells
%     %% cell radius Initialization
%     % Cell radius about normal distribution
     Cradius(cells, 1) = randgaussrad(Cell_radius, Cell_std);
% 
% 
%     %% Position Initialization
%     %Polar coordinates to confine system position initialization
    r(cells, 1) = R_boundary*sqrt(rand());
    theta(cells, 1) = 2*pi*rand();
% 
%     % Returns in cartesian components
    x(cells, 1) = lbox/2 + r(cells,1) * cos(theta(cells, 1));
    y(cells, 1) = lbox/2 + r(cells,1) * sin(theta(cells, 1));
    %% Directionality / velocity  initialization
    %Velocity vector angle initialization
    vel_ang(cells, 1) = 2*pi*rand();
    % Cell speed about normal distribution
    %speed(cells, 1) = randgaussrad(vels_med, vels_std);
    speed(cells,1) = vels_med;
    %Split componentwise
    vx(cells,1)= speed(cells, 1)*cos(vel_ang(cells,1));
    vy(cells,1)= speed(cells,1)*sin(vel_ang(cells,1));
end % end for loop
end % end function