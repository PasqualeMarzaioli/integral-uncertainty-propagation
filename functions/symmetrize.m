% symmetrize.m — Enforce symmetry of a covariance matrix.
% Replaces P with (P + P^T)/2 to remove roundoff-level antisymmetric
% terms before Cholesky or eigendecomposition.
% Author: Pasquale Marzaioli

function covariance = symmetrize(covariance)
% Remove roundoff-level antisymmetric terms before decompositions.
covariance = (covariance + covariance.') / 2;
end
