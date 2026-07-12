% transformCovarianceSeries.m — Rotate covariance histories into RTH.
% For each nominal apsis mean, builds the RTH frame and applies the
% block-diagonal transform T = blkdiag(R,R) to LinCov, UT, and Monte Carlo
% covariances: P_rth = T P_eci T^T.
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
