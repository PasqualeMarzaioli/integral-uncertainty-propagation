% unscentedSigmaPoints.m — Construct scaled unscented-transform sigma points.
% Returns the deterministic points and their separate mean and covariance weights.
% Author: Pasquale Marzaioli

function [points, meanWeights, covarianceWeights] = unscentedSigmaPoints( ...
    meanState, covariance, alpha, beta, kappa)
% Construct the scaled sigma set and its two weight vectors.
dimension = numel(meanState);
lambda = alpha^2 * (dimension + kappa) - dimension;
scale = dimension + lambda;
squareRoot = chol(scale * covariance, "lower");
points = [meanState, meanState + squareRoot, meanState - squareRoot];
meanWeights = [lambda / scale; repmat(1 / (2 * scale), 2 * dimension, 1)];
covarianceWeights = meanWeights;
covarianceWeights(1) = covarianceWeights(1) + 1 - alpha^2 + beta;
end
