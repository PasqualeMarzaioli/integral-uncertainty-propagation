% variationalDynamics.m — Propagate the state and state-transition matrix.
% Uses the analytic two-body gravity gradient to evolve the variational equations.
% Author: Pasquale Marzaioli

function derivative = variationalDynamics(time, augmentedState, mu)
% Propagate the nominal state and the 6-by-6 state-transition matrix.
state = augmentedState(1:6);
position = state(1:3);
radius = norm(position);
gravityGradient = mu * (3 * (position * position.') / radius^5 ...
    - eye(3) / radius^3);
dynamicsJacobian = [zeros(3), eye(3); gravityGradient, zeros(3)];
stateTransitionMatrix = reshape(augmentedState(7:end), 6, 6);
derivative = [twoBodyDynamics(time, state, mu); ...
    reshape(dynamicsJacobian * stateTransitionMatrix, 36, 1)];
end
