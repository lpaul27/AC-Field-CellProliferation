% Function for electric field forces

% Inputs: cell position(x,y), E field componenents, E-field grid
% Outputs: component Force exerted on cell(i)

function [EF_x, EF_y, Epressure] = Electric_Force( Cradius, x, y, u, v, X, Y)

% From initialization parameter in EF_Grid_Init function
global NumCells                                                             %#ok<GVMIS> 
EF_x = zeros(NumCells, 1);
EF_y = zeros(NumCells, 1);

for i = 1:NumCells
    % Loops over each cell and find the relevant grid points to avoid
    % calulation of unnecessary points
    
    % Filters the range to only use points which are relevant to the field
    % interaction
    coord_sub = find(X>(x(i)-Cradius(i)) & X<(x(i)+Cradius(i)) & Y>(y(i)-Cradius(i)) & Y<(y(i)+Cradius(i)));              
    
    % Force is calculated based on the range of important values to the
    % cell
    EF_x(i) = sum(u(coord_sub)) ./ length(u(coord_sub));
    EF_y(i) = sum(v(coord_sub)) ./ length(v(coord_sub));
end
% calculate "electric field pressure"
Epressure = sqrt(EF_x.^2 + EF_y.^2) ./ (pi*Cradius);

end


