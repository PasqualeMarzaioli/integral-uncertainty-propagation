% variationalDynamics.m — Nominal state plus 6x6 STM ODE.
% Propagates the augmented vector [x; vec(Phi)] where Phi is the state
% transition matrix. The gravity-gradient Jacobian for two-body motion is
%   A = [0 I; G 0],  G = mu (3 r r^T / r^5 - I / r^3),
% and Phi_dot = A Phi with Phi(0) = I. This is the core of linear
% covariance (LinCov) propagation.
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
