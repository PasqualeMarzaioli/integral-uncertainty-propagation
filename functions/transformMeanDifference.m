% transformMeanDifference.m — Rotate Cartesian mean differences into RTH.
% Applies the same block-diagonal RTH transform used for covariances to a
% sequence of six-component mean differences (e.g. UT minus LinCov).
% Author: Pasquale Marzaioli

function transformed = transformMeanDifference(difference, frames)
% Rotate a sequence of six-component Cartesian differences into RTH.
transformed = zeros(size(difference));
for index = 1:size(difference, 2)
    transform = blkdiag(frames(:, :, index), frames(:, :, index));
    transformed(:, index) = transform * difference(:, index);
end
end
