% plotMonteCarloStability.m — Plot cumulative Monte Carlo mean stability.
% Compares final-pericentre and final-apocentre component histories against the stable index.
% Author: Pasquale Marzaioli

function plotMonteCarloStability(sampleCount, pericentreMean, apocentreMean, ...
    minimumStableSampleCount, figureName, fileName, saveFlag, folder)
% Show the cumulative sample mean at the final pericentre and apocentre.
figure("Name", figureName, "Color", "w");
layout = tiledlayout(6, 2, "TileSpacing", "compact");
stateLabels = ["r_1 / km", "r_2 / km", "r_3 / km", ...
    "v_1 / km/s", "v_2 / km/s", "v_3 / km/s"];
for component = 1:6
    nexttile;
    plot(sampleCount, pericentreMean(component, :), "k");
    xline(minimumStableSampleCount, "--", "HandleVisibility", "off");
    grid on;
    ylabel(stateLabels(component));
    if component == 1
        title("Final pericentre");
    end
    nexttile;
    plot(sampleCount, apocentreMean(component, :), "k");
    xline(minimumStableSampleCount, "--", "HandleVisibility", "off");
    grid on;
    if component == 1
        title("Final apocentre");
    end
end
xlabel(layout, "Number of Monte Carlo samples");
saveNamedFigure(gcf, fileName, saveFlag, folder);
end
