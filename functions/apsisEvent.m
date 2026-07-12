% apsisEvent.m — ODE event detecting apoapsis/periapsis.
% Locates stationary radius through roots of r · v = 0 (radial velocity
% zero). isTerminal=0 and direction=0 so every apsis is reported without
% stopping the integration.
% Author: Pasquale Marzaioli

function [value, isTerminal, direction] = apsisEvent(~, state, ~)
% Detect stationary radius through roots of r dot v.
value = dot(state(1:3), state(4:6));
isTerminal = 0;
direction = 0;
end
