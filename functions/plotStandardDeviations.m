% plotStandardDeviations.m — LinCov vs UT one-sigma in RTH.
% Compares diagonal standard deviations of LinCov and unscented
% covariances for every RTH position and velocity component.
% Author: Pasquale Marzaioli

function plotStandardDeviations(linear, unscented, apsisLabel, figureName, fileName, saveFlag, folder)
% Compare LinCov and UT one-sigma histories for every RTH state component.
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
    grid on;
    ylabel(sprintf("\\sigma_{r,%s} / km", componentLabels(component)));
    if component == 1
        legend("LinCov", "UT", "Location", "best");
    end
    nexttile;
    plot(apsisIndices, sqrt(reshape(linear(component + 3, component + 3, :), [], 1)), ...
        "-o", "LineWidth", 1.5);
    hold on;
    plot(apsisIndices, sqrt(reshape(unscented(component + 3, component + 3, :), [], 1)), ...
        "--s", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\sigma_{v,%s} / km/s", componentLabels(component)));
end
xlabel(apsisLabel);
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
