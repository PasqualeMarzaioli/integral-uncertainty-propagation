% plotCovarianceEllipse.m — Draw a scaled two-dimensional covariance contour.
% Converts the covariance eigensystem into a sampled ellipse around the requested center.
% Author: Pasquale Marzaioli

function plotCovarianceEllipse(covariance, center, sigmaScale, color, lineStyle)
% Plot a centered two-dimensional covariance ellipse from its eigensystem.
if nargin < 5
    lineStyle = "-";
end
[eigenvectors, eigenvalues] = eig(symmetrize(covariance));
angles = linspace(0, 2 * pi, 240);
ellipse = sigmaScale * eigenvectors * sqrt(max(eigenvalues, 0)) ...
    * [cos(angles); sin(angles)] + center;
plot(ellipse(1, :), ellipse(2, :), "Color", color, ...
    "LineStyle", lineStyle, "LineWidth", 1.5);
end
