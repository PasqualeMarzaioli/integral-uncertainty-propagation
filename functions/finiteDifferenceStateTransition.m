% finiteDifferenceStateTransition.m — Numerical STM via central differences.
% Independently verifies the analytically integrated Phi by perturbing
% each initial-state component and re-propagating to targetTime.
% Mixed steps [1e-3 km; 1e-7 km/s] balance position and velocity scales.
% Author: Pasquale Marzaioli

function matrix = finiteDifferenceStateTransition(initialState, targetTime, mu, options)
% Independently verify the propagated STM using scaled central differences.
steps = [1e-3; 1e-3; 1e-3; 1e-7; 1e-7; 1e-7];
matrix = zeros(6);
for component = 1:6
    perturbation = zeros(6, 1);
    perturbation(component) = steps(component);
    upperState = propagateState(initialState + perturbation, targetTime, mu, options);
    lowerState = propagateState(initialState - perturbation, targetTime, mu, options);
    matrix(:, component) = (upperState(:, end) - lowerState(:, end)) ...
        / (2 * steps(component));
end
end
