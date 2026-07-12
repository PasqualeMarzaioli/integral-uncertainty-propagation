% assertCovarianceSeries.m — Assert symmetry and PSD of a covariance series.
% Checks every epoch of a 6x6xN history for Frobenius symmetry and
% positive semidefiniteness (allowing tiny negative eigenvalues).
% Author: Pasquale Marzaioli

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
