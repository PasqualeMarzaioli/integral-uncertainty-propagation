% transformCovarianceSeries.m — Rotate covariance histories into local RTH frames.
% Uses each nominal reference state to rotate all three propagation methods consistently.
% Author: Pasquale Marzaioli

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
