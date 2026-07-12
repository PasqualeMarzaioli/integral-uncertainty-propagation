% plotMaximumCovarianceAxes.m — Max 3-sigma axes across LinCov/UT/MC.
% For pericentre and apocentre series, plots the largest eigenvalue-based
% 3-sigma position and velocity axes and their differences from Monte Carlo.
% Author: Pasquale Marzaioli

function plotMaximumCovarianceAxes(linearPericentre, unscentedPericentre, monteCarloPericentre, ...
    linearApocentre, unscentedApocentre, monteCarloApocentre, ...
    figureName, fileName, saveFlag, folder)
% Compare three-sigma major axes and their differences from Monte Carlo.
methodsPericentre = {linearPericentre, unscentedPericentre, monteCarloPericentre};
methodsApocentre = {linearApocentre, unscentedApocentre, monteCarloApocentre};
labels = ["LinCov", "UT", "MC"];
figure("Name", figureName, "Color", "w");
tiledlayout(4, 2, "TileSpacing", "compact");
for column = 1:2
    if column == 1
        methodSet = methodsPericentre;
        apsisLabel = "Pericentre";
    else
        methodSet = methodsApocentre;
        apsisLabel = "Apocentre";
    end
    count = size(methodSet{1}, 3);
    apsisIndices = 1:count;
    positionAxes = zeros(3, count);
    velocityAxes = zeros(3, count);
    for method = 1:3
        for index = 1:count
            positionAxes(method, index) = 3 * sqrt(max(eig( ...
                symmetrize(methodSet{method}(1:3, 1:3, index)))));
            velocityAxes(method, index) = 3 * sqrt(max(eig( ...
                symmetrize(methodSet{method}(4:6, 4:6, index)))));
        end
    end
    nexttile(column);
    plot(apsisIndices, positionAxes, "LineWidth", 1.4);
    grid on;
    ylabel("3\sigma_{r,max} / km");
    if column == 1
        legend(labels, "Location", "best");
    end
    nexttile(column + 2);
    plot(apsisIndices, positionAxes(1:2, :) - positionAxes(3, :), "LineWidth", 1.4);
    grid on;
    ylabel("Difference / km");
    nexttile(column + 4);
    plot(apsisIndices, velocityAxes, "LineWidth", 1.4);
    grid on;
    ylabel("3\sigma_{v,max} / km/s");
    nexttile(column + 6);
    plot(apsisIndices, velocityAxes(1:2, :) - velocityAxes(3, :), "LineWidth", 1.4);
    grid on;
    ylabel("Difference / km/s");
    xlabel(apsisLabel);
end
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
