% saveNamedFigure.m — Export a figure under a stable publication name.
% Creates the output directory on demand and removes interactive axis toolbars.
% Author: Pasquale Marzaioli

function saveNamedFigure(figureHandle, fileName, saveFlag, folder)
% Export only when requested; descriptive names remain stable in the code.
if ~saveFlag
    return;
end
if ~isfolder(folder)
    mkdir(folder);
end
set(findall(figureHandle, "Type", "axes"), "Toolbar", []);
exportgraphics(figureHandle, fullfile(folder, fileName + ".png"), ...
    "Resolution", 300);
end
