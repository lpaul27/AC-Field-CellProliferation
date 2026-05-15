% Function for polarity Case Handling
function [polarity, vx, vy] = Polarization(polarity, Ex, Ey, vels_ang, time)

global E0 tau vels_med runTime dt

if(Ex >= E0)
    E_ratio = 1;
end
if(Ex < E0 )
    E_ratio = Ex/E0;
end
field_ang = atan2(Ey,Ex);
polarity(:,1) = polarity(:,1) + (E_ratio - polarity(:,1)) * dt* tau / (runTime); 
vx(:,1) = vels_med .* (1+E_ratio).*((1-polarity(:,1)).* cos(vels_ang(:,1)) - polarity(:,1) .* cos(field_ang));
vy(:,1) = vels_med .* (1+E_ratio).*((1-polarity(:,1)).* sin(vels_ang(:,1)) - polarity(:,1) .* sin(field_ang));

end