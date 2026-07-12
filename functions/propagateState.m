% propagateState.m — Integrate one realization to multiple target times.
% Uses ode113 with twoBodyDynamics from t=0 to every requested epoch;
% drops the initial sample so the output columns match targetTimes.
% Author: Pasquale Marzaioli

function states = propagateState(initialState, targetTimes, mu, options)
% Integrate one independent realization to every requested epoch.
[~, history] = ode113(@(time, state) twoBodyDynamics(time, state, mu), ...
    [0; targetTimes], initialState, options);
states = history(2:end, :).';
end
