% symmetrize.m — Remove roundoff-level matrix antisymmetry.
% Averages a matrix with its transpose before covariance decompositions.
% Author: Pasquale Marzaioli

function covariance = symmetrize(covariance)
% Remove roundoff-level antisymmetric terms before decompositions.
covariance = (covariance + covariance.') / 2;
end
