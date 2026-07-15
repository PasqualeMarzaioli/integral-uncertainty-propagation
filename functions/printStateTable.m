% printStateTable.m — Print apsis epochs and Cartesian states as aligned text.
% Preserves enough solver digits for deterministic comparison with the results document.
% Author: Pasquale Marzaioli

function printStateTable(epochs, states, label)
% Print apsis states with the same significant-digit style used in the report.
headers = ["Apsis", "EpochUTC", "rx_km", "ry_km", "rz_km", ...
    "vx_km_s", "vy_km_s", "vz_km_s"];
rows = [repmat(string(label), numel(epochs), 1), string(epochs(:)), ...
    compose("%#.15g", states.')];
widths = max(strlength([headers; rows]), [], 1) + 2;
for column = 1:numel(headers)
    fprintf("%-*s", widths(column), char(headers(column)));
end
fprintf("\n");
for row = 1:size(rows, 1)
    for column = 1:numel(headers)
        fprintf("%-*s", widths(column), char(rows(row, column)));
    end
    fprintf("\n");
end
fprintf("\n");
end
