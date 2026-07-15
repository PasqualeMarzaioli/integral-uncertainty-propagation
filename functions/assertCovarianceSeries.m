% assertCovarianceSeries.m — Validate a propagated covariance history.
% Checks finite values, symmetry, variances, and dimensionless positive semidefiniteness.
% Author: Pasquale Marzaioli

function assertCovarianceSeries(series, methodName)
% Check finite values, symmetry, and PSD after removing mixed state units.
for index = 1:size(series, 3)
    covariance = series(:, :, index);
    assert(all(isfinite(covariance), "all"), ...
        "%s covariance is non-finite at epoch %d.", methodName, index);
    variances = diag(covariance);
    assert(all(variances >= -1e-12 * max(abs(variances), realmin)), ...
        "%s covariance has a negative variance at epoch %d.", methodName, index);
    standardDeviations = sqrt(max(variances, realmin));
    normalizedCovariance = covariance ./ ...
        (standardDeviations * standardDeviations.');
    symmetryError = norm(normalizedCovariance - normalizedCovariance.', "fro");
    assert(symmetryError <= 1e-10, ...
        "%s covariance lost symmetry at epoch %d.", methodName, index);
    minimumEigenvalue = min(eig(symmetrize(normalizedCovariance)));
    assert(minimumEigenvalue >= -1e-10, ...
        "%s covariance is not positive semidefinite at epoch %d.", methodName, index);
end
end
