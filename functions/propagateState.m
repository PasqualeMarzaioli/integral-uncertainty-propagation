% propagateState.m — Propagate one state to requested output epochs.
% Integrates the two-body dynamics with the caller-provided solver options.
% Author: Pasquale Marzaioli

function states = propagateState(initialState, targetTimes, mu, options)
% Integrate one independent realization to every requested epoch.
[~, history] = ode113(@(time, state) twoBodyDynamics(time, state, mu), ...
    [0; targetTimes], initialState, options);
states = history(2:end, :).';
end
