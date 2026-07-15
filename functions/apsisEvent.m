% apsisEvent.m — Detect orbital apsides during numerical integration.
% Locates roots of radial motion through the scalar product of position and velocity.
% Author: Pasquale Marzaioli

function [value, isTerminal, direction] = apsisEvent(~, state, ~)
% Detect stationary radius through roots of r dot v.
value = dot(state(1:3), state(4:6));
isTerminal = 0;
direction = 0;
end
