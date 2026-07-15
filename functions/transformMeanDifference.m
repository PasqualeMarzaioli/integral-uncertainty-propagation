% transformMeanDifference.m — Rotate Cartesian mean differences into RTH.
% Applies the nominal local frame independently at each propagated epoch.
% Author: Pasquale Marzaioli

function transformed = transformMeanDifference(difference, frames)
% Rotate a sequence of six-component Cartesian differences into RTH.
transformed = zeros(size(difference));
for index = 1:size(difference, 2)
    transform = blkdiag(frames(:, :, index), frames(:, :, index));
    transformed(:, index) = transform * difference(:, index);
end
end
