% plotMeanDifference.m — Plot UT-minus-LinCov mean histories in RTH coordinates.
% Shows every position and velocity component plus their Euclidean norms.
% Author: Pasquale Marzaioli

function plotMeanDifference(difference, apsisLabel, figureName, fileName, saveFlag, folder)
% Plot UT-minus-LinCov component and norm histories in local RTH coordinates.
figure("Name", figureName, "Color", "w");
tiledlayout(4, 2, "TileSpacing", "compact");
componentLabels = ["r", "t", "h"];
apsisIndices = 1:size(difference, 2);
for component = 1:3
    nexttile;
    plot(apsisIndices, difference(component, :), "-ok", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\Delta r_{%s} / km", componentLabels(component)));
    nexttile;
    plot(apsisIndices, difference(component + 3, :), "-ok", "LineWidth", 1.5);
    grid on;
    ylabel(sprintf("\\Delta v_{%s} / km/s", componentLabels(component)));
end
nexttile;
plot(apsisIndices, vecnorm(difference(1:3, :), 2, 1), "-ok", "LineWidth", 1.5);
grid on;
ylabel("||\Delta r|| / km");
xlabel(apsisLabel);
nexttile;
plot(apsisIndices, vecnorm(difference(4:6, :), 2, 1), "-ok", "LineWidth", 1.5);
grid on;
ylabel("||\Delta v|| / km/s");
xlabel(apsisLabel);
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
