% twoBodyDynamics.m — Evaluate the two-body state derivative.
% Computes point-mass Earth acceleration in the benchmark inertial frame.
% Author: Pasquale Marzaioli

function derivative = twoBodyDynamics(~, state, mu)
% Evaluate point-mass Earth dynamics in the benchmark inertial frame.
position = state(1:3);
radius = norm(position);
derivative = [state(4:6); -mu * position / radius^3];
end
