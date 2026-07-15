% plotStandardDeviations.m — Compare RTH component standard deviations.
% Draws LinCov, UT, and Monte Carlo histories with distinct line and marker styles.
% Author: Pasquale Marzaioli

function plotStandardDeviations(linear, unscented, monteCarlo, ...
    apsisLabel, figureName, fileName, saveFlag, folder)
% Compare all three one-sigma histories for every RTH state component.
figure("Name", figureName, "Color", "w");
tiledlayout(3, 2, "TileSpacing", "compact");
componentLabels = ["r", "t", "h"];
apsisIndices = 1:size(linear, 3);
for component = 1:3
    nexttile;
    plot(apsisIndices, sqrt(reshape(linear(component, component, :), [], 1)), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, sqrt(reshape(unscented(component, component, :), [], 1)), ...
        "--s", "LineWidth", 1.5);
    plot(apsisIndices, sqrt(reshape(monteCarlo(component, component, :), [], 1)), ...
        ":d", "Color", "k", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\sigma_{r,%s} / km", componentLabels(component)));
    if component == 1
        legend("LinCov", "UT", "MC", "Location", "best");
    end
    nexttile;
    plot(apsisIndices, sqrt(reshape(linear(component + 3, component + 3, :), [], 1)), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, sqrt(reshape(unscented(component + 3, component + 3, :), [], 1)), ...
        "--s", "LineWidth", 1.5);
    plot(apsisIndices, sqrt(reshape( ...
        monteCarlo(component + 3, component + 3, :), [], 1)), ...
        ":d", "Color", "k", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\sigma_{v,%s} / km/s", componentLabels(component)));
end
xlabel(apsisLabel);
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
