% Function for interaction forces
% Adhesion and repulsion

% Inputs: Position of all cells, radius of all cells, total number of
% cells, angle of travel wrt norm
% Outputs: componentwise cell-cell interaction force, average angle of
% neighbors within reaction radius

function [Fx, Fy, neibAngAvg, Pressure] = Interaction_Forces(x, y, Cradius, vel_ang, exempt)

% Constants in function
global k NumCells adh neighborWeight %#ok<GVMIS>

%% Raw Computations
% Define meshgrid to quantify overlap, radius, death
[XX,YY] = meshgrid(x, y);
[RadGrid] = meshgrid(Cradius);
[exemptGrid] = meshgrid(exempt);

% Distance between cells
sepx = XX.' - XX;
sepy = YY - YY.';
dist_btw_cell = sqrt(sepx.*sepx + sepy.*sepy);
sep_angle = atan2(sepy,sepx);

% Define grid of cell radius
sum_cell_radii = RadGrid.'+RadGrid;

% computes overlap unfiltered by true overlap
overlap_raw = sum_cell_radii - dist_btw_cell;
%% Repulsive forces [vectorized]

% Define grid based on filtered overlap
logicalGrid = overlap_raw > 0  & overlap_raw < sum_cell_radii;

% filters dead cells
logicalGrid = logicalGrid .* exemptGrid';
trueOverlap = overlap_raw.*logicalGrid;
anglesep = sep_angle.*logicalGrid;

% Calculate force of repulsion
repulsionx =  -k*trueOverlap.*(sepx./(dist_btw_cell+eye(NumCells)));
repulsiony = -k*trueOverlap.*(sepy./(dist_btw_cell+eye(NumCells)));
Frx = sum((repulsionx),1)';
Fry = sum(repulsiony,1)';

%% Adhesion Forces [vectorized]
% Apply grid to adhesion

% center to center distance
cellRij = dist_btw_cell.*logicalGrid;

% radius of cell 'i'
cellRi = logicalGrid.* RadGrid;

% radius of cell 'j'
cellRj = logicalGrid .* RadGrid';

% simplify calculations by terms
term1 = (2.*cellRi.*cellRij).^2;                        
term2 = (cellRi.^2 - cellRj.^2 + cellRij.^2).^2; 

% **avoids sqrt(-1) in odd cases**
diff = abs(term1 - term2); 

%'lij' calculation (vertical ovcerlap length)
cell_vert_overlap = sqrt(diff)./cellRij;

% NaN avoidance
cell_vert_overlap(isnan(cell_vert_overlap)) = 0; 

% Calculate force of adhesion based on model
adhesionx = adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* cos(anglesep);
adhesiony = adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* sin(anglesep);
Fax =  (sum(adhesionx))';
Fay = (sum(adhesiony))';
%% Calculate Net quantities

% Forces
Fx = Fax + Frx;
Fy = Fay + Fry;

% Pressure
Force_gridx = (adhesionx) + repulsionx;
Force_gridy = (adhesiony) + repulsiony;

% 1D pressure
Pressure = sqrt(Force_gridx.^2 + Force_gridy.^2) ./ cell_vert_overlap;
Pressure(isnan(Pressure)) = 0;

Pressure = (sum(Pressure))';
%% Neighboring Group Allignment
% find which cells are within interaction radius

[angleGrid] = meshgrid(vel_ang);
% dynamic allignment radius vector
alignment_radius = 20 * Cradius;

% Sort by grid
index_grid = (dist_btw_cell <= alignment_radius & dist_btw_cell > 0);

% weight all cell ij interactions (i=j unweighted)
index_grid = index_grid * (neighborWeight/ (1+neighborWeight)) + eye(size(index_grid));

% convert to cartesian to find average
angleGridX = cos(angleGrid);
angleGridY = sin(angleGrid);
angGridXT = angleGridX';
angGridYT = angleGridY';

% sort to only necessary values
InteractionRadGridX = (index_grid .* angGridXT);
InteractionRadGridY = (index_grid .* angGridYT);

% Take componentwise arctan to get average angle within a radius
AngSumX = sum(InteractionRadGridX);
AngSumY = sum(InteractionRadGridY);
neibAngAvg = (atan2(AngSumY, AngSumX))';
end % end function