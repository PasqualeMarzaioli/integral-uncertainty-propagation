% integral_uncertainty_propagation.m — Propagate INTEGRAL state uncertainty.
% The script detects apsides and independently compares LinCov, UT, and Monte Carlo propagation.
% Author: Pasquale Marzaioli

clear;
close all;
clc;

%% Configuration and initial estimate

muEarth = 398600.4418; % km^3/s^2
initialEpoch = datetime(2025, 10, 26, 1, 10, 57, 769, ...
    "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ss.SSS");
initialMean = [7896.37847595; -3152.23575783; 149961.25963557; ...
              -0.6035374932; 0.0159618300; 0.0321154658];
initialCovariance = [6.7e-3, -2.1e-3,  7.6e-4, 0, 0, 0; ...
                    -2.1e-3, 9.7e-3,  5.3e-4, 0, 0, 0; ...
                     7.6e-4, 5.3e-4,  2.1e-3, 0, 0, 0; ...
                     0,       0,       0,       4.7e-7, 0, 0; ...
                     0,       0,       0,       0, 3.1e-7, 0; ...
                     0,       0,       0,       0, 0, 1.6e-7];

numberOfApsides = 5;
numberOfConvergenceSamples = 2000;
convergenceToleranceSigma = 0.05;
randomSeed = 42;
saveFigures = true;
figureDirectory = "plots";

integratorOptions = odeset("RelTol", 1e-12, "AbsTol", 1e-13);

%% Detect the next five pericentres and apocentres

initialRadius = norm(initialMean(1:3));
initialSpeedSquared = dot(initialMean(4:6), initialMean(4:6));
semimajorAxis = 1 / (2 / initialRadius - initialSpeedSquared / muEarth);
orbitalPeriod = 2 * pi * sqrt(semimajorAxis^3 / muEarth);

eventOptions = odeset(integratorOptions, "Events", ...
    @(time, state) apsisEvent(time, state, muEarth));
propagationDuration = (numberOfApsides + 0.1) * orbitalPeriod;
[~, ~, eventTimes, eventStates] = ode113( ...
    @(time, state) twoBodyDynamics(time, state, muEarth), ...
    [0, propagationDuration], initialMean, eventOptions);

% Ignore a possible root at the supplied initial apocentre, then classify
% each later root through the time derivative of r dot v at the event.
validEvents = eventTimes > 1;
eventTimes = eventTimes(validEvents);
eventStates = eventStates(validEvents, :);
radialTurningDerivative = sum(eventStates(:, 4:6).^2, 2) ...
    - muEarth ./ vecnorm(eventStates(:, 1:3), 2, 2);
pericentreMask = radialTurningDerivative > 0;
apocentreMask = radialTurningDerivative < 0;

pericentreTimes = eventTimes(pericentreMask);
pericentreStates = eventStates(pericentreMask, :).';
apocentreTimes = eventTimes(apocentreMask);
apocentreStates = eventStates(apocentreMask, :).';

assert(numel(pericentreTimes) >= numberOfApsides, ...
    "Fewer than five pericentre events were detected.");
assert(numel(apocentreTimes) >= numberOfApsides, ...
    "Fewer than five apocentre events were detected.");
pericentreTimes = pericentreTimes(1:numberOfApsides);
pericentreStates = pericentreStates(:, 1:numberOfApsides);
apocentreTimes = apocentreTimes(1:numberOfApsides);
apocentreStates = apocentreStates(:, 1:numberOfApsides);

pericentreEpochs = initialEpoch + seconds(pericentreTimes);
apocentreEpochs = initialEpoch + seconds(apocentreTimes);
fprintf("Detected apsides for five revolutions (Keplerian period %.3f s).\n", orbitalPeriod);
printStateTable(pericentreEpochs, pericentreStates, "Pericentre");
printStateTable(apocentreEpochs, apocentreStates, "Apocentre");

%% Linearized covariance propagation using the STM (LinCov)

[targetTimes, targetOrder] = sort([pericentreTimes; apocentreTimes]);
pericentreTargetIndices = find(targetOrder <= numberOfApsides);
apocentreTargetIndices = find(targetOrder > numberOfApsides);

initialAugmentedState = [initialMean; reshape(eye(6), 36, 1)];
[~, augmentedHistory] = ode113( ...
    @(time, state) variationalDynamics(time, state, muEarth), ...
    [0; targetTimes], initialAugmentedState, integratorOptions);
augmentedHistory = augmentedHistory(2:end, :).';

linearMeans = augmentedHistory(1:6, :);
linearCovariances = zeros(6, 6, numel(targetTimes));
for epochIndex = 1:numel(targetTimes)
    stateTransitionMatrix = reshape(augmentedHistory(7:end, epochIndex), 6, 6);
    linearCovariances(:, :, epochIndex) = symmetrize( ...
        stateTransitionMatrix * initialCovariance * stateTransitionMatrix.');
end

%% Unscented-transform propagation

alpha = 0.1;
beta = 2;
kappa = 0;
[sigmaPoints, meanWeights, covarianceWeights] = unscentedSigmaPoints( ...
    initialMean, initialCovariance, alpha, beta, kappa);

numberOfSigmaPoints = size(sigmaPoints, 2);
propagatedSigmaPoints = zeros(6, numel(targetTimes), numberOfSigmaPoints);
for sigmaIndex = 1:numberOfSigmaPoints
    propagatedSigmaPoints(:, :, sigmaIndex) = propagateState( ...
        sigmaPoints(:, sigmaIndex), targetTimes, muEarth, integratorOptions);
end

unscentedMeans = zeros(6, numel(targetTimes));
unscentedCovariances = zeros(6, 6, numel(targetTimes));
for epochIndex = 1:numel(targetTimes)
    epochSigmaPoints = reshape(propagatedSigmaPoints(:, epochIndex, :), 6, []);
    unscentedMeans(:, epochIndex) = epochSigmaPoints * meanWeights;
    centeredSigmaPoints = epochSigmaPoints - unscentedMeans(:, epochIndex);
    unscentedCovariances(:, :, epochIndex) = symmetrize( ...
        (centeredSigmaPoints .* covarianceWeights.') * centeredSigmaPoints.');
end

%% Monte Carlo propagation and convergence

rng(randomSeed, "twister");
initialSquareRoot = chol(initialCovariance, "lower");
initialSamples = initialMean + initialSquareRoot * ...
    randn(6, numberOfConvergenceSamples);
propagatedSamples = zeros(6, numel(targetTimes), numberOfConvergenceSamples);

fprintf("Propagating %d independent Monte Carlo samples...\n", ...
    numberOfConvergenceSamples);
progressStep = max(1, floor(numberOfConvergenceSamples / 10));
for sampleIndex = 1:numberOfConvergenceSamples
    propagatedSamples(:, :, sampleIndex) = propagateState( ...
        initialSamples(:, sampleIndex), targetTimes, muEarth, integratorOptions);
    if mod(sampleIndex, progressStep) == 0
        fprintf("  Monte Carlo progress: %d/%d\n", ...
            sampleIndex, numberOfConvergenceSamples);
    end
end

analysisSamples = propagatedSamples;
monteCarloMeans = zeros(6, numel(targetTimes));
monteCarloCovariances = zeros(6, 6, numel(targetTimes));
for epochIndex = 1:numel(targetTimes)
    epochSamples = reshape(analysisSamples(:, epochIndex, :), 6, []);
    monteCarloMeans(:, epochIndex) = mean(epochSamples, 2);
    monteCarloCovariances(:, :, epochIndex) = cov(epochSamples.');
end

finalPericentreSamples = reshape( ...
    propagatedSamples(:, pericentreTargetIndices(end), :), 6, []);
finalApocentreSamples = reshape( ...
    propagatedSamples(:, apocentreTargetIndices(end), :), 6, []);
sampleCount = 1:numberOfConvergenceSamples;
pericentreCumulativeMean = cumsum(finalPericentreSamples, 2) ./ sampleCount;
apocentreCumulativeMean = cumsum(finalApocentreSamples, 2) ./ sampleCount;

% Define convergence as remaining within a fixed fraction of each component's
% final sample standard deviation for every subsequent cumulative estimate.
pericentreScale = max(std(finalPericentreSamples, 0, 2), eps);
apocentreScale = max(std(finalApocentreSamples, 0, 2), eps);
pericentreNormalizedError = max(abs(pericentreCumulativeMean ...
    - pericentreCumulativeMean(:, end)) ./ pericentreScale, [], 1);
apocentreNormalizedError = max(abs(apocentreCumulativeMean ...
    - apocentreCumulativeMean(:, end)) ./ apocentreScale, [], 1);
remainingNormalizedError = fliplr(cummax(fliplr(max( ...
    pericentreNormalizedError, apocentreNormalizedError))));
minimumStableSampleCount = find( ...
    remainingNormalizedError <= convergenceToleranceSigma, 1);
assert(~isempty(minimumStableSampleCount), ...
    "Monte Carlo convergence could not be established.");
fprintf("Monte Carlo means remain within %.3f sigma after N = %d samples.\n", ...
    convergenceToleranceSigma, minimumStableSampleCount);

%% Transform results to the local RTH frames

linearPericentreMeans = linearMeans(:, pericentreTargetIndices);
linearApocentreMeans = linearMeans(:, apocentreTargetIndices);
unscentedPericentreMeans = unscentedMeans(:, pericentreTargetIndices);
unscentedApocentreMeans = unscentedMeans(:, apocentreTargetIndices);

linearPericentreCovariances = linearCovariances(:, :, pericentreTargetIndices);
linearApocentreCovariances = linearCovariances(:, :, apocentreTargetIndices);
unscentedPericentreCovariances = unscentedCovariances(:, :, pericentreTargetIndices);
unscentedApocentreCovariances = unscentedCovariances(:, :, apocentreTargetIndices);
monteCarloPericentreCovariances = monteCarloCovariances(:, :, pericentreTargetIndices);
monteCarloApocentreCovariances = monteCarloCovariances(:, :, apocentreTargetIndices);

[pericentreFrames, linearPericentreCovariancesRth, ...
    unscentedPericentreCovariancesRth, monteCarloPericentreCovariancesRth] = ...
    transformCovarianceSeries(linearPericentreMeans, ...
    linearPericentreCovariances, unscentedPericentreCovariances, ...
    monteCarloPericentreCovariances);
[apocentreFrames, linearApocentreCovariancesRth, ...
    unscentedApocentreCovariancesRth, monteCarloApocentreCovariancesRth] = ...
    transformCovarianceSeries(linearApocentreMeans, ...
    linearApocentreCovariances, unscentedApocentreCovariances, ...
    monteCarloApocentreCovariances);

pericentreMeanDifferenceRth = transformMeanDifference( ...
    unscentedPericentreMeans - linearPericentreMeans, pericentreFrames);
apocentreMeanDifferenceRth = transformMeanDifference( ...
    unscentedApocentreMeans - linearApocentreMeans, apocentreFrames);

%% Reproduce the uncertainty-propagation figures

plotMeanDifference(pericentreMeanDifferenceRth, "Pericentre", ...
    "UT minus LinCov mean at pericentres", ...
    "pericentre_mean_difference", saveFigures, figureDirectory);
plotMeanDifference(apocentreMeanDifferenceRth, "Apocentre", ...
    "UT minus LinCov mean at apocentres", ...
    "apocentre_mean_difference", saveFigures, figureDirectory);

plotStandardDeviations(linearPericentreCovariancesRth, ...
    unscentedPericentreCovariancesRth, "Pericentre", ...
    "RTH standard deviations at pericentres", ...
    "pericentre_rth_standard_deviations", saveFigures, figureDirectory);
plotStandardDeviations(linearApocentreCovariancesRth, ...
    unscentedApocentreCovariancesRth, "Apocentre", ...
    "RTH standard deviations at apocentres", ...
    "apocentre_rth_standard_deviations", saveFigures, figureDirectory);

plotCrossCovariances(linearPericentreCovariancesRth, ...
    unscentedPericentreCovariancesRth, "Pericentre", ...
    "RTH cross-covariances at pericentres", ...
    "pericentre_rth_cross_covariances", saveFigures, figureDirectory);
plotCrossCovariances(linearApocentreCovariancesRth, ...
    unscentedApocentreCovariancesRth, "Apocentre", ...
    "RTH cross-covariances at apocentres", ...
    "apocentre_rth_cross_covariances", saveFigures, figureDirectory);

plotMonteCarloConvergence(sampleCount, pericentreCumulativeMean, ...
    apocentreCumulativeMean, minimumStableSampleCount, ...
    "Monte Carlo mean convergence", ...
    "monte_carlo_mean_convergence", saveFigures, figureDirectory);

plotMonteCarloEllipses(analysisSamples(:, pericentreTargetIndices, :), ...
    linearPericentreMeans, unscentedPericentreMeans, ...
    monteCarloMeans(:, pericentreTargetIndices), ...
    linearPericentreCovariances, unscentedPericentreCovariances, ...
    monteCarloPericentreCovariances, pericentreFrames, "Pericentre", ...
    "Monte Carlo distributions at pericentres", ...
    "pericentre_monte_carlo_distributions", saveFigures, figureDirectory);
plotMonteCarloEllipses(analysisSamples(:, apocentreTargetIndices, :), ...
    linearApocentreMeans, unscentedApocentreMeans, ...
    monteCarloMeans(:, apocentreTargetIndices), ...
    linearApocentreCovariances, unscentedApocentreCovariances, ...
    monteCarloApocentreCovariances, apocentreFrames, "Apocentre", ...
    "Monte Carlo distributions at apocentres", ...
    "apocentre_monte_carlo_distributions", saveFigures, figureDirectory);

plotMaximumCovarianceAxes(linearPericentreCovariances, ...
    unscentedPericentreCovariances, monteCarloPericentreCovariances, ...
    linearApocentreCovariances, unscentedApocentreCovariances, ...
    monteCarloApocentreCovariances, "Maximum covariance axes", ...
    "maximum_covariance_axes", saveFigures, figureDirectory);

%% Physical and numerical verification

assert(max(abs(diff(pericentreTimes) - orbitalPeriod)) < 0.1, ...
    "Pericentre spacing is inconsistent with the Keplerian period.");
assert(max(abs(diff(apocentreTimes) - orbitalPeriod)) < 0.1, ...
    "Apocentre spacing is inconsistent with the Keplerian period.");

allNominalStates = [pericentreStates, apocentreStates];
specificEnergy = vecnorm(allNominalStates(4:6, :), 2, 1).^2 / 2 ...
    - muEarth ./ vecnorm(allNominalStates(1:3, :), 2, 1);
referenceEnergy = initialSpeedSquared / 2 - muEarth / initialRadius;
assert(max(abs((specificEnergy - referenceEnergy) / referenceEnergy)) < 1e-10, ...
    "Specific orbital energy was not conserved to the requested tolerance.");

firstStateTransitionMatrix = reshape(augmentedHistory(7:end, 1), 6, 6);
finiteDifferenceMatrix = finiteDifferenceStateTransition( ...
    initialMean, targetTimes(1), muEarth, integratorOptions);
relativeStmError = norm(firstStateTransitionMatrix - finiteDifferenceMatrix, "fro") ...
    / norm(finiteDifferenceMatrix, "fro");
assert(relativeStmError < 1e-5, ...
    "The variational STM disagrees with a central finite-difference check.");

assertCovarianceSeries(linearCovariances, "LinCov");
assertCovarianceSeries(unscentedCovariances, "UT");
assertCovarianceSeries(monteCarloCovariances, "Monte Carlo");
fprintf("Uncertainty-propagation verification passed: apsis timing, energy, symmetry, and PSD checks.\n");

%% Local functions

function derivative = twoBodyDynamics(~, state, mu)
% Evaluate point-mass Earth dynamics in ECI J2000.
position = state(1:3);
radius = norm(position);
derivative = [state(4:6); -mu * position / radius^3];
end

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

function [value, isTerminal, direction] = apsisEvent(~, state, ~)
% Detect stationary radius through roots of r dot v.
value = dot(state(1:3), state(4:6));
isTerminal = 0;
direction = 0;
end

function states = propagateState(initialState, targetTimes, mu, options)
% Integrate one independent realization to every requested epoch.
[~, history] = ode113(@(time, state) twoBodyDynamics(time, state, mu), ...
    [0; targetTimes], initialState, options);
states = history(2:end, :).';
end

function [points, meanWeights, covarianceWeights] = unscentedSigmaPoints( ...
    meanState, covariance, alpha, beta, kappa)
% Construct the scaled Julier-Uhlmann sigma set and its two weight vectors.
dimension = numel(meanState);
lambda = alpha^2 * (dimension + kappa) - dimension;
scale = dimension + lambda;
squareRoot = chol(scale * covariance, "lower");
points = [meanState, meanState + squareRoot, meanState - squareRoot];
meanWeights = [lambda / scale; repmat(1 / (2 * scale), 2 * dimension, 1)];
covarianceWeights = meanWeights;
covarianceWeights(1) = covarianceWeights(1) + 1 - alpha^2 + beta;
end

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

function frame = rthFrame(state)
% Build row basis vectors that map ECI components into radial-transverse-normal.
radial = state(1:3) / norm(state(1:3));
normal = cross(state(1:3), state(4:6));
normal = normal / norm(normal);
transverse = cross(normal, radial);
frame = [radial.'; transverse.'; normal.'];
end

function [frames, linearRth, unscentedRth, monteCarloRth] = ...
    transformCovarianceSeries(referenceMeans, linear, unscented, monteCarlo)
% Rotate three covariance histories into each nominal apsis RTH frame.
count = size(referenceMeans, 2);
frames = zeros(3, 3, count);
linearRth = zeros(6, 6, count);
unscentedRth = zeros(6, 6, count);
monteCarloRth = zeros(6, 6, count);
for index = 1:count
    frames(:, :, index) = rthFrame(referenceMeans(:, index));
    transform = blkdiag(frames(:, :, index), frames(:, :, index));
    linearRth(:, :, index) = transform * linear(:, :, index) * transform.';
    unscentedRth(:, :, index) = transform * unscented(:, :, index) * transform.';
    monteCarloRth(:, :, index) = transform * monteCarlo(:, :, index) * transform.';
end
end

function transformed = transformMeanDifference(difference, frames)
% Rotate a sequence of six-component Cartesian differences into RTH.
transformed = zeros(size(difference));
for index = 1:size(difference, 2)
    transform = blkdiag(frames(:, :, index), frames(:, :, index));
    transformed(:, index) = transform * difference(:, index);
end
end

function printStateTable(epochs, states, label)
% Print apsis states with the same significant-digit style used in the report.
headers = ["Apsis", "EpochUTC", "rx_km", "ry_km", "rz_km", ...
    "vx_km_s", "vy_km_s", "vz_km_s"];
rows = [repmat(string(label), numel(epochs), 1), string(epochs(:)), ...
    compose("%#.15g", states.')];
widths = max(strlength([headers; rows]), [], 1) + 2;
for column = 1:numel(headers)
    fprintf("%-*s", widths(column), char(headers(column)));
end
fprintf("\n");
for row = 1:size(rows, 1)
    for column = 1:numel(headers)
        fprintf("%-*s", widths(column), char(rows(row, column)));
    end
    fprintf("\n");
end
fprintf("\n");
end

function plotMeanDifference(difference, apsisLabel, figureName, fileName, saveFlag, folder)
% Plot UT-minus-LinCov component and norm histories in local RTH coordinates.
figure("Name", figureName, "Color", "w");
tiledlayout(4, 2, "TileSpacing", "compact");
componentLabels = ["r", "t", "h"];
apsisIndices = 1:size(difference, 2);
for component = 1:3
    nexttile;
    plot(apsisIndices, difference(component, :), "-ok", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\Delta r_{%s} / km", componentLabels(component)));
    nexttile;
    plot(apsisIndices, difference(component + 3, :), "-ok", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\Delta v_{%s} / km/s", componentLabels(component)));
end
nexttile;
plot(apsisIndices, vecnorm(difference(1:3, :), 2, 1), "-ok", "LineWidth", 1.5);
grid on;
ylabel("||\Delta r|| / km");
xlabel(apsisLabel);
nexttile;
plot(apsisIndices, vecnorm(difference(4:6, :), 2, 1), "-ok", "LineWidth", 1.5);
grid on;
ylabel("||\Delta v|| / km/s");
xlabel(apsisLabel);
saveNamedFigure(gcf, fileName, saveFlag, folder);
end

function plotStandardDeviations(linear, unscented, apsisLabel, figureName, fileName, saveFlag, folder)
% Compare LinCov and UT one-sigma histories for every RTH state component.
figure("Name", figureName, "Color", "w");
tiledlayout(3, 2, "TileSpacing", "compact");
componentLabels = ["r", "t", "h"];
apsisIndices = 1:size(linear, 3);
for component = 1:3
    nexttile;
    plot(apsisIndices, sqrt(reshape(linear(component, component, :), [], 1)), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, sqrt(reshape(unscented(component, component, :), [], 1)), ...
        "--s", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\sigma_{r,%s} / km", componentLabels(component)));
    if component == 1
        legend("LinCov", "UT", "Location", "best");
    end
    nexttile;
    plot(apsisIndices, sqrt(reshape(linear(component + 3, component + 3, :), [], 1)), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, sqrt(reshape(unscented(component + 3, component + 3, :), [], 1)), ...
        "--s", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\sigma_{v,%s} / km/s", componentLabels(component)));
end
xlabel(apsisLabel);
saveNamedFigure(gcf, fileName, saveFlag, folder);
end

function plotCrossCovariances(linear, unscented, apsisLabel, figureName, fileName, saveFlag, folder)
% Compare the three distinct within-position and within-velocity covariances.
figure("Name", figureName, "Color", "w");
tiledlayout(3, 2, "TileSpacing", "compact");
pairs = [1, 2; 1, 3; 2, 3];
pairLabels = ["rt", "rh", "th"];
apsisIndices = 1:size(linear, 3);
for pairIndex = 1:3
    row = pairs(pairIndex, 1);
    column = pairs(pairIndex, 2);
    nexttile;
    plot(apsisIndices, reshape(linear(row, column, :), [], 1), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, reshape(unscented(row, column, :), [], 1), ...
        "--s", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("P_{r,%s} / km^2", pairLabels(pairIndex)));
    if pairIndex == 1
        legend("LinCov", "UT", "Location", "best");
    end
    nexttile;
    plot(apsisIndices, reshape(linear(row + 3, column + 3, :), [], 1), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, reshape(unscented(row + 3, column + 3, :), [], 1), ...
        "--s", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("P_{v,%s} / km^2/s^2", pairLabels(pairIndex)));
end
xlabel(apsisLabel);
saveNamedFigure(gcf, fileName, saveFlag, folder);
end

function plotMonteCarloConvergence(sampleCount, pericentreMean, apocentreMean, ...
    minimumStableSampleCount, figureName, fileName, saveFlag, folder)
% Show the cumulative sample mean at the final pericentre and apocentre.
figure("Name", figureName, "Color", "w");
tiledlayout(6, 2, "TileSpacing", "compact");
stateLabels = ["r_1 / km", "r_2 / km", "r_3 / km", ...
    "v_1 / km/s", "v_2 / km/s", "v_3 / km/s"];
for component = 1:6
    nexttile;
    plot(sampleCount, pericentreMean(component, :), "k");
    xline(minimumStableSampleCount, "--", "N_{min}", "HandleVisibility", "off");
    grid on;
    ylabel(stateLabels(component));
    nexttile;
    plot(sampleCount, apocentreMean(component, :), "k");
    xline(minimumStableSampleCount, "--", "N_{min}", "HandleVisibility", "off");
    grid on;
end
xlabel("Number of Monte Carlo samples");
saveNamedFigure(gcf, fileName, saveFlag, folder);
end

function plotMonteCarloEllipses(samples, linearMean, unscentedMean, monteCarloMean, ...
    linearCovariance, unscentedCovariance, monteCarloCovariance, frames, ...
    apsisLabel, figureName, fileName, saveFlag, folder)
% Project samples, means, and three-sigma covariance ellipses into the R-T plane.
figure("Name", figureName, "Color", "w");
count = size(samples, 2);
tiledlayout(1, count, "TileSpacing", "compact");
for index = 1:count
    nexttile;
    frame = frames(:, :, index);
    transform = blkdiag(frame, frame);
    epochSamples = reshape(samples(:, index, :), 6, []);
    sampleDisplacements = frame * (epochSamples(1:3, :) - linearMean(1:3, index));
    plot(sampleDisplacements(1, :), sampleDisplacements(2, :), ".", ...
        "Color", [0.65, 0.65, 0.65], "MarkerSize", 3);
    hold on;
    linearOffset = frame * (linearMean(1:3, index) - linearMean(1:3, index));
    unscentedOffset = frame * (unscentedMean(1:3, index) - linearMean(1:3, index));
    monteCarloOffset = frame * (monteCarloMean(1:3, index) - linearMean(1:3, index));
    plotCovarianceEllipse(transform(1:2, :) * linearCovariance(:, :, index) ...
        * transform(1:2, :).', linearOffset(1:2), 3, [0, 0.4470, 0.7410]);
    plotCovarianceEllipse(transform(1:2, :) * unscentedCovariance(:, :, index) ...
        * transform(1:2, :).', unscentedOffset(1:2), 3, [0.8500, 0.3250, 0.0980]);
    plotCovarianceEllipse(transform(1:2, :) * monteCarloCovariance(:, :, index) ...
        * transform(1:2, :).', monteCarloOffset(1:2), 3, [0, 0, 0], "--");
    plot(linearOffset(1), linearOffset(2), "+", "LineWidth", 1.5, "MarkerSize", 8);
    plot(unscentedOffset(1), unscentedOffset(2), "+", "LineWidth", 1.5, "MarkerSize", 8);
    plot(monteCarloOffset(1), monteCarloOffset(2), "+k", "LineWidth", 1.5, "MarkerSize", 8);
    axis equal;
    grid on;
    horizontalLimits = xlim;
    xticks(linspace(horizontalLimits(1), horizontalLimits(2), 3));
    xlabel("r_r / km");
    if index == 1
        ylabel("r_t / km");
        legend("Samples", "LinCov", "UT", "MC", "Location", "best");
    end
    title(sprintf("%s %d", apsisLabel, index));
end
saveNamedFigure(gcf, fileName, saveFlag, folder);
end

function plotCovarianceEllipse(covariance, center, sigmaScale, color, lineStyle)
% Plot a centered two-dimensional covariance ellipse from its eigensystem.
if nargin < 5
    lineStyle = "-";
end
[eigenvectors, eigenvalues] = eig(symmetrize(covariance));
angles = linspace(0, 2 * pi, 240);
ellipse = sigmaScale * eigenvectors * sqrt(max(eigenvalues, 0)) ...
    * [cos(angles); sin(angles)] + center;
plot(ellipse(1, :), ellipse(2, :), "Color", color, ...
    "LineStyle", lineStyle, "LineWidth", 1.5);
end

function plotMaximumCovarianceAxes(linearPericentre, unscentedPericentre, monteCarloPericentre, ...
    linearApocentre, unscentedApocentre, monteCarloApocentre, ...
    figureName, fileName, saveFlag, folder)
% Compare three-sigma major axes and their differences from Monte Carlo.
methodsPericentre = {linearPericentre, unscentedPericentre, monteCarloPericentre};
methodsApocentre = {linearApocentre, unscentedApocentre, monteCarloApocentre};
labels = ["LinCov", "UT", "MC"];
figure("Name", figureName, "Color", "w");
tiledlayout(4, 2, "TileSpacing", "compact");
for column = 1:2
    if column == 1
        methodSet = methodsPericentre;
        apsisLabel = "Pericentre";
    else
        methodSet = methodsApocentre;
        apsisLabel = "Apocentre";
    end
    count = size(methodSet{1}, 3);
    apsisIndices = 1:count;
    positionAxes = zeros(3, count);
    velocityAxes = zeros(3, count);
    for method = 1:3
        for index = 1:count
            positionAxes(method, index) = 3 * sqrt(max(eig( ...
                symmetrize(methodSet{method}(1:3, 1:3, index)))));
            velocityAxes(method, index) = 3 * sqrt(max(eig( ...
                symmetrize(methodSet{method}(4:6, 4:6, index)))));
        end
    end
    nexttile(column);
    plot(apsisIndices, positionAxes, "LineWidth", 1.4);
    grid on;
    ylabel("3\sigma_{r,max} / km");
    if column == 1
        legend(labels, "Location", "best");
    end
    nexttile(column + 2);
    plot(apsisIndices, positionAxes(1:2, :) - positionAxes(3, :), "LineWidth", 1.4);
    grid on;
    ylabel("Difference / km");
    nexttile(column + 4);
    plot(apsisIndices, velocityAxes, "LineWidth", 1.4);
    grid on;
    ylabel("3\sigma_{v,max} / km/s");
    nexttile(column + 6);
    plot(apsisIndices, velocityAxes(1:2, :) - velocityAxes(3, :), "LineWidth", 1.4);
    grid on;
    ylabel("Difference / km/s");
    xlabel(apsisLabel);
end
saveNamedFigure(gcf, fileName, saveFlag, folder);
end

function covariance = symmetrize(covariance)
% Remove roundoff-level antisymmetric terms before decompositions.
covariance = (covariance + covariance.') / 2;
end

function assertCovarianceSeries(series, methodName)
% Check symmetry and positive semidefiniteness at every propagated epoch.
for index = 1:size(series, 3)
    covariance = series(:, :, index);
    symmetryError = norm(covariance - covariance.', "fro");
    scale = max(1, norm(covariance, "fro"));
    assert(symmetryError <= 1e-10 * scale, ...
        "%s covariance lost symmetry at epoch %d.", methodName, index);
    minimumEigenvalue = min(eig(symmetrize(covariance)));
    assert(minimumEigenvalue >= -1e-10 * scale, ...
        "%s covariance is not positive semidefinite at epoch %d.", methodName, index);
end
end

function saveNamedFigure(figureHandle, fileName, saveFlag, folder)
% Export only when requested; descriptive names remain stable in the code.
if ~saveFlag
    return;
end
if ~isfolder(folder)
    mkdir(folder);
end
exportgraphics(figureHandle, fullfile(folder, fileName + ".png"), ...
    "Resolution", 300);
end
