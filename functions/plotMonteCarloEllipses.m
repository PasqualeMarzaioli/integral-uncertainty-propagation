% plotMonteCarloEllipses.m — MC samples and 3-sigma ellipses in the R-T plane.
% For each apsis, projects samples relative to the LinCov mean into the
% local radial-transverse plane and overlays LinCov, UT, and MC ellipses.
% Author: Pasquale Marzaioli

function plotMonteCarloEllipses(samples, linearMean, unscentedMean, monteCarloMean, ...
    linearCovariance, unscentedCovariance, monteCarloCovariance, frames, ...
    apsisLabel, figureName, fileName, saveFlag, folder)
% Project samples, means, and three-sigma covariance ellipses into the R-T plane.
figure("Name", figureName, "Color", "w");
count = size(samples, 2);
tiledlayout(1, count, "TileSpacing", "compact");
for index = 1:count
    nexttile;
    frame = frames(:, :, index);
    transform = blkdiag(frame, frame);
    epochSamples = reshape(samples(:, index, :), 6, []);
    sampleDisplacements = frame * (epochSamples(1:3, :) - linearMean(1:3, index));
    plot(sampleDisplacements(1, :), sampleDisplacements(2, :), ".", ...
        "Color", [0.65, 0.65, 0.65], "MarkerSize", 3);
    hold on;
    linearOffset = frame * (linearMean(1:3, index) - linearMean(1:3, index));
    unscentedOffset = frame * (unscentedMean(1:3, index) - linearMean(1:3, index));
    monteCarloOffset = frame * (monteCarloMean(1:3, index) - linearMean(1:3, index));
    plotCovarianceEllipse(transform(1:2, :) * linearCovariance(:, :, index) ...
        * transform(1:2, :).', linearOffset(1:2), 3, [0, 0.4470, 0.7410]);
    plotCovarianceEllipse(transform(1:2, :) * unscentedCovariance(:, :, index) ...
        * transform(1:2, :).', unscentedOffset(1:2), 3, [0.8500, 0.3250, 0.0980]);
    plotCovarianceEllipse(transform(1:2, :) * monteCarloCovariance(:, :, index) ...
        * transform(1:2, :).', monteCarloOffset(1:2), 3, [0, 0, 0], "--");
    plot(linearOffset(1), linearOffset(2), "+", "LineWidth", 1.5, "MarkerSize", 8);
    plot(unscentedOffset(1), unscentedOffset(2), "+", "LineWidth", 1.5, "MarkerSize", 8);
    plot(monteCarloOffset(1), monteCarloOffset(2), "+k", "LineWidth", 1.5, "MarkerSize", 8);
    axis equal;
    grid on;
    horizontalLimits = xlim;
    xticks(linspace(horizontalLimits(1), horizontalLimits(2), 3));
    xlabel("r_r / km");
    if index == 1
        ylabel("r_t / km");
        legend("Samples", "LinCov", "UT", "MC", "Location", "best");
    end
    title(sprintf("%s %d", apsisLabel, index));
end
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
