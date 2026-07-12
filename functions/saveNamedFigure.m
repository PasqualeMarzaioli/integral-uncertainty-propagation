% saveNamedFigure.m — Conditionally export a figure to PNG.
% When saveFlag is true, writes figureHandle into folder as fileName.png
% at 300 dpi. Shared by all three navigation scripts so plot export
% behaviour stays consistent.
% Author: Pasquale Marzaioli

function saveNamedFigure(figureHandle, fileName, saveFlag, folder)
% Export only when requested while keeping descriptive names stable.
if ~saveFlag
    return;
end
if ~isfolder(folder)
    mkdir(folder);
end
exportgraphics(figureHandle, fullfile(folder, string(fileName) + ".png"), ...
    "Resolution", 300);
end
