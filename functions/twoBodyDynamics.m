% twoBodyDynamics.m — Point-mass Earth two-body ODE right-hand side.
% Autonomous Keplerian dynamics in ECI J2000:
%   r_dot = v,  v_dot = -mu r / r^3.
% Used by Monte Carlo and as the nominal part of variationalDynamics.
% Author: Pasquale Marzaioli

function derivative = twoBodyDynamics(~, state, mu)
% Evaluate point-mass Earth dynamics in ECI J2000.
position = state(1:3);
radius = norm(position);
derivative = [state(4:6); -mu * position / radius^3];
end
