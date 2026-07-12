% plotCrossCovariances.m — LinCov vs UT RTH cross-covariances.
% Plots the three distinct within-position and within-velocity off-diagonal
% terms (rt, rh, th) for both methods versus apsis index.
% Author: Pasquale Marzaioli

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
