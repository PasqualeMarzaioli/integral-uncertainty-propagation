% plotMonteCarloConvergence.m — Cumulative MC mean at final apsides.
% 6-by-2 panels of running sample means for the final pericentre and
% apocentre states, with a vertical marker at the minimum stable sample count.
% Author: Pasquale Marzaioli

function plotMonteCarloConvergence(sampleCount, pericentreMean, apocentreMean, ...
    minimumStableSampleCount, figureName, fileName, saveFlag, folder)
% Show the cumulative sample mean at the final pericentre and apocentre.
figure("Name", figureName, "Color", "w");
tiledlayout(6, 2, "TileSpacing", "compact");
stateLabels = ["r_1 / km", "r_2 / km", "r_3 / km", ...
    "v_1 / km/s", "v_2 / km/s", "v_3 / km/s"];
for component = 1:6
    nexttile;
    plot(sampleCount, pericentreMean(component, :), "k");
    xline(minimumStableSampleCount, "--", "N_{min}", "HandleVisibility", "off");
    grid on;
    ylabel(stateLabels(component));
    nexttile;
    plot(sampleCount, apocentreMean(component, :), "k");
    xline(minimumStableSampleCount, "--", "N_{min}", "HandleVisibility", "off");
    grid on;
end
xlabel("Number of Monte Carlo samples");
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
